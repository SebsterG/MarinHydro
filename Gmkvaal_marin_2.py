import numpy as np 
from timeit import default_timer as timer
from math import pi
from matplotlib.pyplot import *


def points(a,b,N):
	# Creating an ellipse/circle
	x = np.zeros(N*2)
	y = np.zeros(N*2)
	theta = 2*pi
	d_theta = theta/N
	for i in range(N*2):
		x[i] = a*np.cos(i*d_theta)
		y[i] = b*np.sin(i*d_theta)
	return x,y

def req_points(b,h,N):
	x = np.zeros(2*N)
	y = np.zeros(2*N)
	x[0:N/4] = b/2.
	x[N/4:N/2] = np.linspace(b/2.,-b/2.,N/4.)
	x[N/2:3*N/4] = -b/2.
	x[3*N/4:N] = np.linspace(-b/2.,b/2.,N/4.)
	y[0:N/4] = np.linspace(-h/2.,h/2.,N/4.)
	y[N/4:N/2] = h/2.
	y[N/2:3*N/4] = np.linspace(h/2.,-h/2.,N/4.)
	y[3*N/4:N] = -h/2.
	x[N:2*N] = x[0:N]
	y[N:2*N] = y[0:N]
	return x,y

x,y = req_points(1,1,100)
plot(x[0:100./2],y[0:100./2])
#show()


def anglefinder(a,b,N):
	# Making a matrix containing the angles between lines from a current point
	# to the joints between the line segments 
	#x,y = points(a,b,N)
	x,y = req_points(a,b,N) 
	matrix = np.zeros([N,N])
	for j in range(N):
		for i in range(N):	
			# Making the vectors and calculating the angle
			v0x = x[i] - 0.5*(x[j]+x[j+1])
			v0y = y[i] - 0.5*(y[j]+y[j+1])
			v1x = x[i+1] - 0.5*(x[j]+x[j+1])
			v1y = y[i+1] - 0.5*(y[j]+y[j+1])
			# Filling out the diagonal, which corresponds to the current point, with -pi 
			if i==j:
				matrix[j,i] = -pi
			else:
				num =  - np.arccos ( (v0x*v1x+v0y*v1y) / 
						( np.sqrt(v0x**2+v0y**2)*np.sqrt(v1x**2+v1y**2)))
				if float('-inf') < float(num) < float('inf'):
					matrix[j,i] = num	
				else:
					matrix[j,i] = 0
	return matrix


def integral(direction,a,b,N,x1,x2,y1,y2):
	#x0,y0 = points(a,b,N)
	x0,y0 = req_points(a,b,N)
	#for i in range(N):
	# Midle of integration segment
	y0c = 0.5*(y0[0:N]+y0[1:N+1])
	x0c = 0.5*(x0[0:N]+x0[1:N+1])
	# Midle of current segment
	x_c = 0.5*(x1+x2)
	y_c	= 0.5*(y1+y2)
	# Distance for current point to endpoints of integration segments 
	rad_A = np.sqrt((x_c-x0[0:N])**2 + (y_c-y0[0:N])**2)
	rad_B = np.sqrt((x_c-x0[1:N+1])**2 + (y_c-y0[1:N+1])**2)
	# Length of integration segment
	ds = np.sqrt((x0[1:N+1]-x0[0:N])**2 + (y0[1:N+1]-y0[0:N])**2)
	if direction == 1:
		# x normal (d_phi_1/dn)
		#n = (- x0c/a**2)/np.sqrt( (x0c**2/a**4)+(y0c**2/b**4) )
		n = (x0[1:N+1]-x0[:N])/ds
		n[np.isnan(n)] = 0
	if direction ==2:
		# y normal (d_phi_2/dn)
		n = (y0[1:N+1]-y0[:N])/ds
		n[np.isnan(n)] = 0
	if direction == 6:
		nx = (x0[1:N+1]-x0[:N])/ds
		nx[np.isnan(nx)] = 0
		ny = (y0[1:N+1]-y0[:N])/ds
		ny[np.isnan(ny)] = 0
		rx = x0[0:N]      
		ry = y0[0:N]
		n = rx*ny - ry*nx
	trick_a = np.nonzero(rad_A)
	trick_b = np.nonzero(rad_B)
 	A = np.log(rad_A[trick_a])
 	B = np.log(rad_B[trick_b])
 	# Trapezoidal rule
	integral = 0.5 * sum(ds[trick_a]*(A+B)*n[trick_a])
	return integral


def list_o_integrals(direction,a,b,N):
	#x,y = points(a,b,N)
	x,y = req_points(a,b,N)
	integral_list = np.zeros(N)
	start = timer()
	for i in range(N):
		integral_list[i] = integral(direction,a,b,N,x[i],x[i+1],y[i],y[i+1])
	end = timer()
	return integral_list


def solver(direction,a,b,N):
	# Solving the matrix equation
	matrix = anglefinder(a,b,N)
	m = matrix
	c = list_o_integrals(direction,a,b,N)
	return np.linalg.solve(m,c)

def reference(r,N):
	d_theta = 2*pi/N
	y = np.zeros(N)
	for i in range(N):
		y[i] = -r*np.cos(i*d_theta)
	return y


def added_mass(direction,a,b,N):
	phi = solver(direction,a,b,N)
	#x0,y0 = points(a,b,N)
	x0,y0 = req_points(a,b,N)
	ds = np.sqrt((x0[1:N+1]-x0[:N])**2 + (y0[1:N+1]-y0[:N])**2)	
	# midpoints
	y0c = 0.5*(y0[:N]+y0[1:N+1])
	x0c = 0.5*(x0[:N]+x0[1:N+1])
	A = phi[:N]
	B = phi[0:N+1]
	print len(A)
	print len(B)
	# Normal x or y component
	if direction == 1:
		#n = -x0c /np.sqrt(x0c**2+y0c**2)
		n = (x0[1:N+1]-x0[:N])/ds
		n[np.isnan(n)] = 0
		ref = pi*b**2
		ref = 1.51*pi*(0.5*b)**2
	if direction == 2:
		#n = -y0c /np.sqrt(x0c**2+y0c**2)
		n = (y0[1:N+1]-y0[:N])/ds
		n[np.isnan(n)] = 0
		ref = 1.51*pi*(0.5*a)**2
	if direction == 6:
		nx = (x0[1:N+1]-x0[:N])/ds
		nx[np.isnan(nx)] = 0
		ny = (y0[1:N+1]-y0[:N])/ds
		ny[np.isnan(ny)] = 0
		rx = x0[0:N]      
		ry = y0[0:N]
		n = rx*ny - ry*nx
		print n
		ref = 0.725*a**4

	integral = 0.5*sum(ds*(A+B)*n)		
	ratio = integral / ref
	print "The added mass coefficient m%d%d of an rectangle with a=%.1f and b=%.1f is %.3f" \
		 % (direction,direction,a,b,integral)
	print "That gives a ratio %.6f with the reference solution %.3f for N=%d" % (ratio,ref,N)


def run_it(direction,a,b,N,Plot,AM):
	start = timer()
	if AM == True:
		added_mass(direction,a,b,N)
	if Plot == True:
		y = solver(direction,a,b,N)
		x = np.linspace(0,360,N)
		amp = abs(min(y)-max(y))
		print "The amplitude is %.3f" % amp 
		plot(x,y)
		title('Potential')
		xlabel('Degrees')
		ylabel('Phi')
		show()
	end = timer()
	print "Total running time: %.3fs" % (end-start)

run_it(direction=6,a=1,b=1,N=1000,Plot=False,AM=True)