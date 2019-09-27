# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 01:20:57 2019

@author: Ankit Raj
"""

# Load necessary modules
import sys
import numpy as np
from scipy import integrate
from scipy.optimize import curve_fit
#-------------------------------------------------------------------
# Define required functions 

#-------------------------------------------------------------------

# Loading data and checks:

distance = np.loadtxt("./data/distance.txt")
potential = np.loadtxt("./data/potential.txt")

correction1 = np.loadtxt("./data/c1.txt")
distance_c1 = np.loadtxt("./data/dc1.txt")

correction2 = np.loadtxt("./data/c2.txt")
distance_c2 = np.loadtxt("./data/dc2.txt")

correction3 = np.loadtxt("./data/c3.txt")
distance_c3 = np.loadtxt("./data/dc3.txt")

#-------------------------------------------------------------------

# Parameters :

# Atomic number
q1 = 1
q1 = 1


# Atomic mass, in atomic units
m1 =  918
m2 = 918


#-------------------------------------------------------------------

