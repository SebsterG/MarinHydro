
from matplotlib.pylab import *
import numpy as np
import math as math


def make_points(num_points,r_a, r_b):
	if r_a == r_b:
		om = 2.0*np.pi*r_a
	else:
		om = np.pi*(3*(r_a+r_b)- np.sqrt(((3*r_a)+r_b)*(r_a+(3*r_b)))) 
	dx = 2*pi/num_points
	x = np.zeros(num_points*3)
	y = np.zeros(num_points*3)
	for i in range(num_points*3):
		x[i] = r_a*np.cos(i*dx)
		y[i] = r_b*np.sin(i*dx)
	return x,y


def make_angle(N,r_a,r_b):
	angle = np.zeros(N-1)
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


	return matrix_
"""	
def make_matrix_A(N,r_a,r_b):
	x,y = make_points(N,r_a,r_b)	
	matrix_ = np.zeros((N,N))
	for i in range(N):
		angle = make_angle(N,r_a,r_b,x[i],x[i+1],y[i],y[i+1])
		for j in range(N):
			if i==j:
				matrix_[i][j] = -np.pi
			else :						
				matrix_[i][j] = -angle[j-1] 

	return matrix_
	"""
#print make_matrix_A(6,1,1)


"""def make_matrix_A(N,r_a,r_b):	
	matrix_ = np.zeros((N,N))
	angle = make_angle(N,r_a,r_b)
	for i in range(N):
		for j in range(N):
			if i==j:
				matrix_[i][j] = -np.pi
				
				index = -j-1
			else :						
				matrix_[i][j] = -angle[j+index] 
	return matrix_
"""
def integral_1(N,r_a,r_b):
	#ds = 2*np.pi*r_a/N
	x,y = make_points(N,r_a,r_b)
	integ = np.zeros(N)
	traps1 = 0
	traps2 = 0
	Matrix_B = np.zeros(N)
	for i in range(N-1):	
		integral_2 = 0
		x0 = x[i]
		y0 = y[i]
		for j in range(N-2):
			rad_1 = np.sqrt((x0-x[j])**2 + (y0-y[j])**2)
			rad_2 = np.sqrt((x0-x[j+1])**2 + (y0-y[j+1])**2)
			if rad_1==0:
				continue
			if rad_2== 0:
				continue
		
			traps1 = (-((x[j]+x[j+1])/2) / (np.sqrt(((x[j]+x[j+1])/2)**2 + ((y[j]+y[j+1])/2)**2))) * np.log(rad_1)	
			traps2 = (-((x[j+1]+x[j+2])/2) / (np.sqrt(((x[j+1]+x[j+2])/2)**2 +((y[j+1]+y[j+2])/2)**2))) * np.log(rad_2)
			ds = np.sqrt((x[j+1]-x[j])**2 + (y[j+1] -y[j])**2)
			integral_2 = integral_2 + (traps1+traps2) * ds * 0.5
		Matrix_B[i] = integral_2
	return Matrix_B
def solve(N,r_a,r_b):
	return np.linalg.solve(make_angle(N,r_a,r_b),integral_1(N,r_a,r_b))

N = 401
r_a = 1
r_b = 1


thet = np.zeros(N)
om = 2.0*pi*r_a
dx = 2*pi/N
for i in range(N):
	thet[i] = i*dx
diff = max(-np.cos(thet) - solve(N,r_a,r_b))
print 100*diff/max(-np.cos(thet))
plot(-np.cos(thet),"g")
legend("exact solution")
plot(solve(N,r_a,r_b), "r")
legend("Distribution of potential over a circle")
show()
