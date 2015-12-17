
import matplotlib.pyplot as plt
import numpy as np
import math as math


def make_points(num_points,r_a, r_b):
	dx = 2*np.pi/num_points
	x = np.zeros(num_points*3)
	y = np.zeros(num_points*3)
	for i in range(num_points*2+1):
		x[i] = r_a*np.cos(i*dx)
		y[i] = r_b*np.sin(i*dx)
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
"""x,y =make_points_1(100,1,1)
plt.figure()
plt.plot(x,y, "-o")
plt.axis([-3,3,-3,3])	
plt.show()"""

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
			if (np.isnan(matrix_[i][k])):
				matrix_[i][k] = 0
	return matrix_


def integral_1(N,r_a,r_b,direction1):
	x,y = make_points(N,r_a,r_b)
	traps1 = 0
	traps2 = 0
	Matrix_B = np.zeros(N)
	ds = np.sqrt((x[1:-1]-x[:-2])**2 + (y[1:-1] - y[:-2])**2)
	for i in range(N):	
		integral_x = 0
		integral_66 = 0
		x0 = 0.5*(x[i]+x[i+1])
		y0 = 0.5*(y[i]+y[i+1])
		#for j in range(N-1):
		rad_1 = np.sqrt((x0-x[:-2])**2 + (y0-y[:-2])**2)
		rad_2 = np.sqrt((x0-x[1:-1])**2 + (y0-y[1:-1])**2)	

		if direction1 == 11:
			print (x[:-2]+x[1:-1])
			traps1x = -((x[:-3]+x[1:-2])/(2*r_a**2)) / np.sqrt( ((x[:-3]+x[1:-2])/(2*r_a**2))**2+\
			((y[:-3]+y[1:-2])/(2*r_b**2))**2  ) * np.log(rad_1[:-1])		
			traps2x = -((x[1:-2]+x[2:-1])/(2*r_a**2)) / np.sqrt(  ((x[1:-2]+x[2:-1])/(2*r_a**2))**2 + ((y[1:-2]+y[2:-1])/(2*r_b**2))**2  )* np.log(rad_2)	
			integral_x = np.sum((traps1x[:] + traps2x[:]) * ds[:])  * 0.5

		else :
			r1a = (x[:-2] + x[1:-1])*0.5 ; r2a =0.5*(y[:-2] + y[1:-1])
			n1a = -((x[:-2]+x[1:-1])/(2*r_a**2))/np.sqrt(  ((x[:-2]+x[1:-1])/(2*r_a**2))**2 + ((y[:-2]+y[1:-1])/(2*r_b**2))**2  )
			n2a = -((y[:-2]+y[1:-1])/(2*r_b**2))/np.sqrt(  ((x[:-2]+x[1:-1])/(2*r_a**2))**2 + ((y[:-2]+y[1:-1])/(2*r_b**2))**2  )
			crosses1 = r1a*n2a - r2a*n1a

			r1b = (x[1:-1] + x[2:])*0.5 ; r2b =0.5*(y[1:-1] + y[2:])
			n1b = -((x[1:-1]+x[2:])/(2*r_a**2))/np.sqrt(  ((x[1:-1]+x[2:])/(2*r_a**2))**2 + ((y[1:-1]+y[2:])/(2*r_b**2))**2  )
			n2b = -((y[1:-1]+y[2:])/(2*r_b**2))/np.sqrt(  ((x[1:-1]+x[2:])/(2*r_a**2))**2 + ((y[1:-1]+y[2:])/(2*r_b**2))**2  )
			crosses2 = r1b*n2b - r2b*n1b

			traps3 = np.log(rad_1) * crosses1
			traps4 = np.log(rad_2) * crosses2
			ds = np.sqrt((x[1:-1]-x[:-2])**2 + (y[1:-1] - y[:-2])**2)
			integral_66 = integral_66 + (traps3 +traps4)* ds   * 0.5
		if direction1  == 11:
			Matrix_B[i] = integral_x
		else :
			Matrix_B[i] = integral_66
	return Matrix_B


def solver1(N,r_a,r_b,direction1):
	return np.linalg.solve(make_angle(N,r_a,r_b),integral_1(N,r_a,r_b,direction1))

def added_mass(N,r_a,r_b):
	phix = solver1(N,r_a,r_b,direction1 =11)
	phi6 = solver1(N,r_a,r_b,direction1 =66)
	x,y = make_points(N,r_a,r_b)
	add_m11 = 0
	add_m22 = 0
	add_m66 = 0
	for i in range(N-2):
		rad_1 = np.sqrt(  ((x[i]+x[i+1])/(2*r_a**2))**2 +  ((y[i] +y[i+1])/(2*r_b**2))**2)
		rad_2 = np.sqrt(  ((x[i+1]+x[i+2])/(2*r_a**2))**2 +  ((y[i+1]+y[i+2])/(2*r_b**2))**2)
		ds = np.sqrt((x[i+1]-x[i])**2 + (y[i+1] - y[i])**2)

		traps1 = phix[i] * -(x[i]+x[i+1])/(2*r_a**2)/rad_1
		traps2 = phix[i+1] * -(x[i+1]+x[i+2])/(2*r_a**2)/rad_2

		r1x = (x[i] + x[i+1])*0.5 ; r2x = 0.5*(y[i] + y[i+1])
		n1x = -(x[i] + x[i+1])/(2*r_a**2)/rad_1
		n2x = -(y[i] + y[i+1])/(2*r_b**2)/rad_1
		crosses1 = r1x*n2x - r2x*n1x

		r1y = (x[i+1] + x[i+2])*0.5 ; r2y = 0.5*(y[i+1] + y[i+2])
		n1y = -(x[i+1] + x[i+2])/(2*r_a**2)/rad_2
		n2y = -(y[i+1] + y[i+2])/(2*r_b**2)/rad_2
		crosses2 = r1y*n2y - r2y*n1y

		traps3 = phi6[i] * crosses1
		traps4 = phi6[i+1] * crosses2

		add_m11 = add_m11 + (traps1 + traps2)*ds * 0.5
		add_m66 = add_m66 + (traps3 + traps4)*ds *0.5
	return add_m11, add_m66




"""
N = 4
r_a = 1
r_b = 1
a11, a66 = added_mass(N,r_a,r_b)
a11_e = 4.754 *r_a**2
a66_e = 0.725 *r_a**2
a11_e = np.pi*r_b**2
#a66_e = (1.0/8.0)* np.pi*(r_a**2-r_b**2)**2
print a11#,a66
print a11_e#, a66_e
#print a11_e, a66_e
"""
N = 1000
r_a = 1
r_b = 1
equal = r_a - r_b

if equal == 0:
	a_m11,a_m66 = added_mass(N,r_a,r_b)
	error1 = ((((r_a**2*np.pi)-float(a_m11))/(r_a**2*np.pi))*100)
	print "m11= %.3f m66= %.f,m11 error in percent: %.2f  %% " %(a_m11,a_m66, error1)
	thet = np.zeros(N)
	om = 2.0*pi*r_a
	dx = 2*pi/N
	for i in range(N):
		thet[i] = i*dx
	solution = solver1(N,r_a,r_b,direction1 = 11)
	exact_1 = -np.cos(thet)
	diff = max(exact_1 - solution)
	plot(exact_1,"g")
	plot(solution, "r")
	plt.legend(['Exact ', 'Numerical'])
	title("Fluid Potential over Circle, running %.d times,  error: %.2f %%  \n m11 = %.2f, m66 = %.f , error in m11 =  %.2f %% " %(N,100*diff/max(exact_1), a_m11 , a_m66, error1)) 
	show()

else :
	a_m11,a_m66 =  added_mass(N,r_a,r_b)
	error1 = (((r_b**2*np.pi)-a_m11)/(r_b**2*np.pi))*100
	exact6 = (np.pi/8)*(r_a**2-r_b**2)**2
	error6 = (((exact6-a_m66 ) / (exact6)))*100
	print "m11: " ,a_m11 ,"error m11 : %.2f  %%, error m66 %.2f %%   running %.d times" %(error1,error6, N)
	plot(solver1(N,r_a,r_b,direction1 = 11), "r")
	title("Distribution of potential over an ellipse running %.d times \n m11: %.2f , error m11 : %.2f %% \n m66 = %.2f , error m66= %.2f %%, r_a = %.d, r_b= %.d  " %(N,a_m11,error1,a_m66,error6,r_a,r_b))
	show()




