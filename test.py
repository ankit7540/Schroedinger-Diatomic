# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 15:53:53 2019

@author: Ankit Raj
"""

import numpy as np

A=np.random.randint(1,50,6)
B=np.random.randint(1,50,(6,6))
print ("A = ",A,"\n \nB = \n", B)

ind = np.argsort(A, axis=0)
print("\nSorting index = \n",ind)

C=A[ind]
print("\nSorted A = \n",C)

D=np.empty_like(B)
for i in range(6):
    D[:,i] = B[:,ind[i]]

print("\nSorted B along the columns, D = \n",D)   


E = B[:,ind] 

print("\nSorted B along the columns, E = \n",E)   

