# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 01:20:57 2019

@author: Ankit Raj
"""

# Load necessary modules
import sys
import numpy as np
from scipy import integrate

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

#-------------------------------------------------------------------
print(start,end)

# Steps

# interpolate  potential and corrections to  the number of  data points

# Generate  the final potential

#  generate the Hmatrix
