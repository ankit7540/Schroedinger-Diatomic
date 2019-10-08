# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 01:20:57 2019

@author: Ankit Raj
"""

# Load necessary modules

#import sys
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
from elements import ELEMENTS



#print  (str(ELEMENTS[8])  )

selection=ELEMENTS[2]
print(selection.number, selection.mass)

m1=selection.isotopes[3]
print(m1)

print( selection.isotopes[4])
print( ELEMENTS[7].isotopes[15])
print( ELEMENTS[8].mass)
print( selection.isotopes[4])

#********************************************************************
# python function to compute the first derivative at the first point
# and the last point of an array. For this computation, the first 4
# points are used for the derivative for the first data point. Similarly
# last 4 (x,y) points are used for the derivative at the last point.

def fd_ends(x, y):
    ''' Parameter:
        x       =       xaxis of the array
        y       =       y data of the array
                (x and y arrays must have atleast 4 elements)
        Returns = first derivative at the first point and the last point
    '''
    if (len(x) < 4 or len(y) < 4):
        print("Error : x and y arrays must have 4 elements")

    subx = np.zeros(4)
    suby = np.zeros(4)

    for i in range(0, 4):
        subx[i] = x[i]
        suby[i] = y[i]

    fd1 = ((subx[1]*subx[3]+subx[2]*subx[3]+subx[1]*subx[2]-2*subx[0]*subx[1] \
          -2*subx[0]*subx[3]-2*subx[0]*subx[2]+3*subx[0]**2)/(-subx[3]+subx[0]) \
          /(-subx[1]+subx[0])/(-subx[2]+subx[0])*suby[0]-(-subx[2]+subx[0])    \
          *(-subx[3]+subx[0])/(subx[1]-subx[3])/(-subx[1]+subx[0])/(subx[1]-subx[2])\
          *suby[1]+(-subx[1]+subx[0])*(-subx[3]+subx[0])/(subx[2]-subx[3])          \
          /(subx[1]-subx[2])/(-subx[2]+subx[0])*suby[2]-(-subx[1]+subx[0])          \
          *(-subx[2]+subx[0])/(subx[2]-subx[3])/(subx[1]-subx[3])/(-subx[3]+subx[0])\
          *suby[3])

    for i in range(0, 4):
        subx[i] = x[int(i-4)]
        suby[i] = y[int(i-4)]
#        print (i, int(i-4))

    fdn = ((subx[1]-subx[3])*(subx[2]-subx[3])/(-subx[3]+subx[0])/(-subx[1] \
           +subx[0])/(-subx[2]+subx[0])*suby[0]-(-subx[3]+subx[0])*(subx[2]   \
           -subx[3])/(subx[1]-subx[3])/(-subx[1]+subx[0])/(subx[1]-subx[2])   \
           *suby[1]+(-subx[3]+subx[0])*(subx[1]-subx[3])/(subx[2]-subx[3])    \
           /(subx[1]-subx[2])/(-subx[2]+subx[0])*suby[2]-(-2*subx[0]*subx[3]  \
           -2*subx[1]*subx[3]-2*subx[2]*subx[3]+subx[0]*subx[1]+subx[0]       \
           *subx[2]+subx[1]*subx[2]+3*subx[3]**2)/(subx[2]-subx[3])/(subx[1]   \
           -subx[3])/(-subx[3]+subx[0])*suby[3])

    return(fd1, fdn)

#************************************************************************
# Define the factorial function
def factorial(n):
    A=np.zeros((5,5))
    C=np.zeros(5)
    C[n]=1

    for(i=0 ; i<5 ; i=i+1)
    	A[][i]=i
    		for (j=0 ; j<5 ; j=j+1)
    			A[j][i]=(i^j) ;
    			A[j][i] = A[j][i] / factorial (j)
    		endfor
    endfor
        return 1
    else:
        return n * factorial(n-1)
#************************************************************************
# Coefs for finite difference derivatives
def coef_r1(n):
    if n == 0:
        return 1
    else:
        return n * factorial(n-1)

#-------------------------------------------------------------------
# Define required functions

#-------------------------------------------------------------------

# Loading data and checks:

distance = np.loadtxt("./data/r_wave670_H2.txt")
potential = np.loadtxt("./data/pes670_H2.txt")

correction1 = np.loadtxt("./data/adiabticE1.txt")
distance_c1 = np.loadtxt("./data/adb_r_distance.txt")

correction2 = np.loadtxt("./data/adiabticE2.txt")
distance_c2 = np.loadtxt("./data/adb_r_distance.txt")

correction3 = np.loadtxt("./data/radiative_57.txt")
distance_c3 = np.loadtxt("./data/r_wave_57_radiative.txt")

#-------------------------------------------------------------------
# Units:
unit_potential = "H"
unit_c1 = "cm-1"
unit_c2 = "cm-1"
unit_c3 = "cm-1"
#-------------------------------------------------------------------

# Parameters :

# Atomic number
q1 = 1
q2 = 1


# Atomic mass, in atomic units
m1 =  918
m2 = 918

#  step  size
step  =  0.005

#  start and end   value of the potential distance curve
start  = distance[0]
end = distance[-1]
rwave=np.arange(start,end,step)
nelements=len(rwave)
#-------------------------------------------------------------------
print(start,end)

# Steps

# interpolate  potential and corrections to  the number of  data points
derv = fd_ends(distance, potential)
cs = CubicSpline(distance, potential,bc_type=((1,derv[0]),(1,derv[1])))
potential_interp = cs(rwave)


derv = fd_ends(distance_c1,  correction1)
cs = CubicSpline(distance_c1,  correction1,bc_type=((1,derv[0]),(1,derv[1])))
adbc1_interp = cs(rwave)


derv = fd_ends(distance_c2,  correction2)
cs = CubicSpline(distance_c2,  correction2,bc_type=((1,derv[0]),(1,derv[1])))
adbc2_interp = cs(rwave)


derv = fd_ends(distance_c3,  correction3)
cs = CubicSpline(distance_c3,  correction3,bc_type=((1,derv[0]),(1,derv[1])))
radc_interp = cs(rwave)


# Generate  the final potential ---------------------------------

#  Generate final potential based on the  unit
if  (unit_potential ==  "cm-1"):
    potential_interp = potential_interp / 219474.631370200000000

if  (unit_c1 ==  "cm-1"):
    adbc1_interp = adbc1_interp / 219474.631370200000000

if  (unit_c2 ==  "cm-1"):
    adbc2_interp = adbc1_interp / 219474.631370200000000

if  (unit_c3 ==  "cm-1"):
    radc_interp = radc_interp / 219474.631370200000000





# final potential
fp = potential_interp+ adbc1_interp+adbc2_interp+radc_interp

#  generate the Hmatrix ------------------------------


H=np.zeros((nelements, nelements), dtype='float')

#for i in range(nelements):
	#print (i )
=======
# coefficients for the derivatives
coefs_d1 = findiff.coefficients(deriv=1, acc=4)
for i in range(5):
    print(i,"\t",coefs_d1['center']['coefficients'][i],"\t",coefs_d1['forward']['coefficients'][i])

print("\n \n")

coefs_d2 = findiff.coefficients(deriv=2, acc=4)
for i in range(5):
    print(i,"\t",coefs_d2['center']['coefficients'][i],"\t",coefs_d2['forward']['coefficients'][i])


# check the coefs with the Igor implementation





H=np.zeros((nelements, nelements), dtype=float)

for i in range(nelements):
	H[i,i]=5




#-------------------------------------------------------------------
#plt.figure(0)
#ax0 = plt.axes()
#plt.title('Potential', fontsize=20)
#plt.plot( rwave, fp ,'r-',  label='potential')
#plt.plot( index,  unexposed_pixel2, 'g-',  label='unexposed  pixel')



#fig, ax = plt.subplots()
# Bilinear interpolation - this will look blurry
#ax1.imshow(H, cmap=cm.RdYlGn )
#plt.show()
