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
import periodictable as pt
#import findiff
#import os
#import dependencies as dep
from dependencies import coefs

#print  (str(ELEMENTS[8])  )



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


#************************************************************************


# Define the factorial function
def factorial(n):
    if n == 0:
        return 1
    else:
        return n * factorial(n-1)

#-------------------------------------------------------------------
# Define required functions

#-------------------------------------------------------------------

# Loading data and checks:

# Born-Oppenheimer potential
distance = np.loadtxt("./data/r_wave670_H2.txt")
potential = np.loadtxt("./data/pes670_H2.txt")

# adiabtaic correction
correction1 = np.loadtxt("./data/adiabticE1.txt")
distance_c1 = np.loadtxt("./data/adb_r_distance.txt")

# adiabtaic correction
correction2 = np.loadtxt("./data/adiabticE2.txt")
distance_c2 = np.loadtxt("./data/adb_r_distance.txt")

# radiative correction
correction3 = np.loadtxt("./data/radiative_57.txt")
distance_c3 = np.loadtxt("./data/r_wave_57_radiative.txt")

# relativistic correction
correction4 = np.loadtxt("./data/relativistic.txt")
distance_c4 = np.loadtxt("./data/r_distance_rel.txt")


#-------------------------------------------------------------------
# Units:
unit_potential = "H"
unit_c1 = "cm-1"
unit_c2 = "cm-1"
unit_c3 = "cm-1"
unit_c4 = "cm-1"
#-------------------------------------------------------------------

# Parameters :

# Atomic number
q1 = 1
q2 = 1



#  step  size
step  =  0.00500

#  start and end   value of the potential distance curve
start  = distance[0]
end = int(distance[-1])
rwave=np.arange(start,end,step)
rwave = np.append(rwave,12.0)
nelements=len(rwave)

hartree_to_wavenumber=2.1947463136320*10**5 
print ("energy:", hartree_to_wavenumber)

print (nelements)
#-------------------------------------------------------------------
#print(start,end)

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


derv = fd_ends(distance_c4,  correction4)
cs = CubicSpline(distance_c4,  correction4, bc_type=((1,derv[0]),(1,derv[1])))
relativistic_interp = cs(rwave)

# Generate  the final potential ---------------------------------

#  Generate final potential based on the  unit
if  (unit_potential ==  "cm-1"):
    potential_interp = potential_interp / hartree_to_wavenumber

if  (unit_c1 ==  "cm-1"):
    adbc1_interp = adbc1_interp / hartree_to_wavenumber

if  (unit_c2 ==  "cm-1"):
    adbc2_interp = adbc1_interp / hartree_to_wavenumber

if  (unit_c3 ==  "cm-1"):
    radc_interp = radc_interp / hartree_to_wavenumber

if  (unit_c4 ==  "cm-1"):
    relativistic_interp = relativistic_interp / hartree_to_wavenumber

# --------------------------------------------------------------


# final potential
fp = potential_interp+ adbc1_interp+adbc2_interp+radc_interp+relativistic_interp


np.savetxt("adiabatic_corr1.txt", adbc1_interp, fmt='%5.12f')
np.savetxt("adiabatic_corr2.txt", adbc2_interp, fmt='%5.12f')
np.savetxt("radiative_corr.txt", radc_interp, fmt='%5.12f')
np.savetxt("relativistiv_corr.txt", relativistic_interp, fmt='%5.12f')


np.savetxt("potential.txt", fp, fmt='%5.12f')
np.savetxt("rwave.txt", rwave, fmt='%5.12f')

#print("interpolated waves generated")
#exit(0)
#  generate the Hmatrix -----------------------------------

def reduced_mass(zA, iMassA,eA, zB, iMassB, eB):
    '''
    zA = atomic number for atom A
    iMassA = atomic mass of specific isotope
    zB = atomic number for atom B
    iMassB = atomic mass of specific isotope
    '''
    mass = 1822.888486209

    sA=pt.elements[zA][iMassA]
    sB=pt.elements[zB][iMassB]

    mA=float(sA.mass)
    mB=float(sB.mass)

    A=mA*mass
    B=mB*mass
    print ("\nA = ",A,"B = ", B)

    #return reduced mass
    return 1/(1/A + 1/B)

#----------------------------------------------------------

def gen_H_matrix(mass,J,rwave,potential, accuracy, step):
    '''
    Generate the Hamiltonian matrix for the radial
    nuclear equation for diatomic molecule
    for fixed accuracy for 5, i.e. 5 point derivative 
    '''
    print ("reduced mass : ", mass)
    H=np.zeros((nelements, nelements), dtype=float)
    ne=int((accuracy-1) / 2)


    for i in range(nelements):

        # assign the diagonal term
        J_term = (J*(J+1)) / (2*mass*(rwave[i])**2 )
        H[i,i]=J_term + potential[i]

        # assign the derivatives

    # first row -------------------------
    D1=coefs.coef_forward (1, accuracy )
    D2=coefs.coef_forward (2, accuracy )

    D1=D1*-1/mass/rwave[0]/step
    D2=D2*-1/(2*mass)/(step**2)

    print (D1)

    for i in range(accuracy):
        H[0][i] = H[0][i] + D1[i] + D2[i]
    #-------------------------------------

    # second row -------------------------
    D1=coefs.coef_forward_asymmetric (1, accuracy, 1 )
    D2=coefs.coef_forward_asymmetric (2, accuracy, 1 )

    D1=D1*-1/mass/rwave[1]/step
    D2=D2*-1/(2*mass)/(step**2)

    for i in range(accuracy):
        H[1][i] = H[1][i] + D1[i] + D2[i]
    #-------------------------------------


    # all symmetric rows -------------------------
    for i in range (2, nelements-2):
        D1=coefs.coef_symmetric_center (1, accuracy )
        D2=coefs.coef_symmetric_center(2, accuracy )

        D1=D1*-1/mass/rwave[i]/step
        D2=D2*-1/(2*mass)/(step**2)
        #print (i, nelements)

        for j in range(accuracy):
            H[i][i-ne+j] = H[i][i-ne+j] + D1[j] + D2[j]
            #print (i,j, i-ne+j, nelements)
    #---------------------------------------------

    # second last row -------------------------
    D1=coefs.coef_backward_asymmetric (1, accuracy, 1 )
    D2=coefs.coef_backward_asymmetric (2, accuracy, 1 )

    D1=D1*-1/mass/rwave[nelements-2]/step
    D2=D2*-1/(2*mass)/(step**2)
    print (nelements-2)

    for i in range(accuracy):
        H[nelements-2][nelements-1-i] = H[nelements-2][nelements-1-i] + D1[i] + D2[i]
    #------------------------------------------

    # last row --------------------------------
    D1=coefs.coef_backward (1, accuracy )
    D2=coefs.coef_backward (2, accuracy )

    D1=D1*-1/mass/rwave[nelements-1]/step
    D2=D2*-1/(2*mass)/(step**2)
    print (nelements-1)

    for i in range(accuracy):
        H[nelements-1][nelements-1-i] = H[nelements-1][nelements-1-i] + D1[i] + D2[i]
    #------------------------------------------




    return H

#-------------------------------------------------------------------
#print(reduced_mass(1,1,1,1,1,1))

#print(reduced_mass(1,1,1,1,2,1))
#print(reduced_mass(1,2,1,1,2,1))
nu = reduced_mass(1,1,1,1,1,1)

H3=gen_H_matrix(nu ,0,rwave,fp,5,step)


#w,v=np.linalg.eig(H3)

#ind = np.argsort(w, axis=0)
#C=w[ind]

#E = w[:,np.argsort(v, axis=0)]

#----------------------------------
H2=H3
H2[0,:]=0
H2[1,:]=0

H2[-1,:]=0
H2[-2,:]=0



#----------------------------------

plt.figure(0)
ax0 = plt.axes()
plt.title('Potential', fontsize=20)
#plt.plot( rwave, E[:,0] ,'r-',  label='potential')
#plt.plot( rwave, E[:,1] ,'r-',  label='potential')
#plt.plot( rwave, E[:,2] ,'r-',  label='potential')



#fig, ax = plt.subplots()
# Bilinear interpolation - this will look blurry
#ax1.imshow(H, cmap=cm.RdYlGn )
#plt.show()
