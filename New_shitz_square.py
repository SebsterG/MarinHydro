from matplotlib.pylab import *
import numpy as np
import math as math

import time
start_time = time.time()
"""
def make_points(num_points,r_a, r_b):
	dx = 2*pi/num_points
	x = np.zeros(num_points*2)
	y = np.zeros(num_points*2)
	for i in range(num_points*2):
		x[i] = r_a*np.cos(i*dx)
		y[i] = r_b*np.sin(i*dx)
	return x,y
"""
"""
def make_points(N,h,b):
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
"""	
def make_points(N,r_a,r_b):
	N = N/4 *4
	Np = N/4
	x = np.zeros(N+1)
	y = np.zeros(N+1)
	d1 = -r_a
	d2 = r_a
	for i in range(Np+1):
		x[i] = d1 +(d2-d1)/2 *(1-np.cos(i*np.pi/Np))
		y[i] = -r_a
		x[i+Np] = r_a
		y[i+Np] = d1 +(d2-d1)/2 *(1-np.cos(i*np.pi/Np))
		x[i+2*Np] = -(d1 +(d2-d1)/2 *(1-np.cos(i*np.pi/Np)))
		y[i+2*Np] = r_a
		x[i+3*Np] = -r_a
		y[i+3*Np] = -(d1 +(d2-d1)/2 *(1-np.cos(i*np.pi/Np)))
	return x ,y

"""
x,y =make_points(100,1,1)
plt.figure()
plt.plot(x,y, "-o")
plt.axis([-3,3,-3,3])	
plt.show()"""

def make_angle(N,r_a,r_b):
	#angle = np.zeros(N-1)
	x,y = make_points(N,r_a,r_b)
	matrix_ = np.zeros((N,N))
	for i in range(N):
		x_first = (x[i]+x[i+1])/2.0
		y_first = (y[i]+y[i+1])/2.0
		for k in range (N):						
			if i == k: 
				matrix_[i][k] = -np.pi
			else:
				xa = x[k] - x_first
				xb = x[k+1] - x_first
				ya = y[k] - y_first
				yb = y[k+1] - y_first
				matrix_[i][k] = -np.arccos((xa*xb + ya*yb)/(np.sqrt(xa**2 + ya**2) * np.sqrt(xb**2 + yb**2)))
			if (np.isnan(matrix_[i][k])):
				matrix_[i][k] = 0
	#print "matrix", matrix_
	return matrix_
#make_angle(100,1,1)
def integral(dirr,N,r_a,r_b,x1,x2,y1,y2):
	x,y = make_points(N,r_a,r_b)
	# middle of this segment
	x_m = (x1+x2)/2.0
	y_m = (y1+y2)/2.0

	# middle of every segment
	# not middle this time
	x_c = (x[1:N+1]-x[:N])
	y_c = (y[1:N+1]-y[:N])
	ds = np.sqrt(x_c**2 + y_c**2)
	#ds = np.sqrt((x[1:N+1]-x[:N])**2 + (y[1:N+1]-y[:N])**2)

	rad_1 = np.sqrt((x_m-x[0:N])**2+(y_m-y[0:N])**2)
	rad_2 = np.sqrt((x_m-x[1:N+1])**2+(y_m-y[1:N+1])**2)

	if dirr == 11:
		n = y_c / ds
		n[np.isnan(n)] = 0
	if dirr == 22:
		n = x_c / ds
		n[np.isnan(n)] = 0
	if dirr == 66:
		ny = x_c / ds
		ny[np.isnan(ny)] = 0
		nx = -y_c / ds
		nx[np.isnan(nx)] = 0
		r_x = (x[1:N+1]+x[:N])/2.0
		r_y = (y[1:N+1]+y[:N])/2.0
		n = ny*r_x - nx*r_y
	"""
	traps1 =np.log(rad_1)
	traps2 = np.log(rad_2)
	"""
	non_1 = np.nonzero(rad_1)
	non_2 = np.nonzero(rad_2)
	traps1 = np.log(rad_1[non_1])
	traps2 = np.log(rad_2[non_2])

	integ = 0.5*sum((traps1+traps2)*ds[non_1]*n[non_1])
	return integ

def integral_matrix(dirr,r_a,r_b,N):
	x,y = make_points(N,r_a,r_b,)
	integral_matrixes = np.zeros(N)
	for i in range(N):
		integral_matrixes[i] = integral(dirr,N,r_a,r_b,x[i],x[i+1],y[i],y[i+1])
	return integral_matrixes
def solver(dirr,N,r_a,r_b):
	return np.linalg.solve(make_angle(N,r_a,r_b),integral_matrix(dirr,r_a,r_b,N))

#plot(solver(11,360,1,1))
#show()

def added_mass(dirr,r_a,r_b,N):
	phi = solver(dirr,N,r_a,r_b)
	x,y = make_points(N,r_a,r_b)
	# middle of every segment
	# not middle this time
	x_c = (x[1:N+1]-x[:N])
	y_c = (y[1:N+1]-y[:N])
	ds = np.sqrt(x_c**2 + y_c**2)
	#ds = np.sqrt((x[1:N+1]-x[:N])**2 + (y[1:N+1]-y[:N])**2)

	
	phi_1 = phi[:N]
	phi_2 = phi[0:N+1]
	if dirr == 11:
		n = y_c / ds
		n[np.isnan(n)] = 0
		#exact = 4.754*r_a**2
		exact = 1.51*np.pi*r_a**2
	if dirr == 22:
		n = x_c / ds
		n[np.isnan(n)] = 0
		#exact = 4.754*r_a**2
		exact = 1.51*np.pi*r_a**2
	if dirr == 66:
		ny = x_c / ds
		ny[np.isnan(ny)] = 0
		nx = -y_c / ds
		nx[np.isnan(nx)] = 0
		r_x = (x[1:N+1]+x[:N])/2.0
		r_y = (y[1:N+1]+y[:N])/2.0
		n = ny*r_x - nx*r_y
		#exact = 0.725*r_a**4
		exact = 0.234*np.pi*r_a**4

	integ = 0.5*sum((phi_1+phi_2)*ds*n)
	print  N," Elements"
	print "Direction :", dirr
	print "Calculated Added mass is: ",integ
	print "Exact added mass:",exact
	print 100-(100*integ/exact),"%"
	print("--- %.2f seconds ---" % (time.time() - start_time))



def exact(r_a,N):
	dtet = 2*pi/N
	value = np.zeros(N)
	for i in range(N):
		value[i] = -r_a*np.cos(dtet*i)
	return value

#print exact(1,10)
#print exact(1,100) - solver(11,100,1,1)

if __name__ == '__main__':
	#added_mass(11,1,1,360)
	
	added_mass(11,1,1,360)
	
	#added_mass(66,1,1,360)

	#added_mass(66,1,1,720)
