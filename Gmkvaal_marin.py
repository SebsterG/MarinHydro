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



def anglefinder(a,b,N):
	# Making a matrix containing the angles between lines from a current point
	# to the joints between the line segments 
	start = timer() 
	#x,y = points(a,b,N)
	x,y = points(a,b,N) 
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
				matrix[j,i] = - np.arccos ( (v0x*v1x+v0y*v1y) / 
				( np.sqrt(v0x**2+v0y**2)*np.sqrt(v1x**2+v1y**2)))

	end	= timer()
	print "The function anglefinder took %.3fs for N=%d" %(end-start,N)
	print matrix
	return matrix
#anglefinder(1,1,100)


def integral(direction,a,b,N,x1,x2,y1,y2):
	#x0,y0 = points(a,b,N)
	x0,y0 = points(a,b,N)
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
		n = (- x0c/a**2)/np.sqrt( (x0c**2/a**4)+(y0c**2/b**4) )
	if direction ==2:
		# y normal (d_phi_2/dn)
		n = - y0c/np.sqrt(x0c**2+y0c**2)
		n = (- y0c/b**2)/np.sqrt( (x0c**2/a**4)+(y0c**2/b**4) )
	if direction == 6:
		nx = (- x0c/a**2)/np.sqrt( (x0c**2/a**4)+(y0c**2/b**4) )
		ny = (- y0c/b**2)/np.sqrt( (x0c**2/a**4)+(y0c**2/b**4) )
		rx = x0c        
		ry = y0c
		n = rx*ny - ry*nx

 	A = np.log(rad_A)
 	B = np.log(rad_B)
 	# Trapezoidal rule
	integral = 0.5 * sum(ds*(A+B)*n)
	return integral


def list_o_integrals(direction,a,b,N):
	x,y = points(a,b,N)
	integral_list = np.zeros(N)
	start = timer()
	for i in range(N):
		integral_list[i] = integral(direction,a,b,N,x[i],x[i+1],y[i],y[i+1])
	end = timer()
	print "list_o_integrals used %.3fs" % (end-start)
	print integral_list
	#integral_list = integral(direction,a,b,N,x[0:N],x[1:N+1],y[0:N],y[1:N+1])	
	return integral_list


def solver(direction,a,b,N):
	# Solving the matrix equation
	#matrix = anglefinder(a,b,N)
	#m = matrix
	#c = list_o_integrals(direction,a,b,N)
	return np.linalg.solve(anglefinder(a,b,N),list_o_integrals(direction,a,b,N))


def reference(r,N):
	d_theta = 2*pi/N
	y = np.zeros(N)
	for i in range(N):
		y[i] = -r*np.cos(i*d_theta)
	return y


def added_mass(direction,a,b,N):
	phi = solver(direction,a,b,N)
	#x0,y0 = points(a,b,N)
	x0,y0 = points(a,b,N)
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
		n = (- x0c/a**2)/np.sqrt( (x0c**2/a**4)+(y0c**2/b**4) )
		ref = pi*b**2
	if direction == 2:
		#n = -y0c /np.sqrt(x0c**2+y0c**2)
		n = (- y0c/b**2)/np.sqrt( (x0c**2/a**4)+(y0c**2/b**4) )
		ref = pi*b**2
	if direction == 6:
		nx = (- x0c/a**2)/np.sqrt( (x0c**2/a**4)+(y0c**2/b**4) )
		ny = (- y0c/b**2)/np.sqrt( (x0c**2/a**4)+(y0c**2/b**4) )
		rx = x0c        
		ry = y0c
		n = rx*ny - ry*nx
		ref = pi*(a**2-b**2)**2/8

	integral = 0.5*sum(ds*(A+B)*n)		
	ratio = integral / ref
	print "The added mass coefficient m%d%d of an ellipse with a=%.1f and b=%.1f is %.3f" \
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

run_it(direction=1,a=1,b=3,N=360,Plot=False,AM=True)

print 1.51*pi*(0.5*2)**2

#print 0.234*pi*(0.5*1)**2