'''
Programs for the coefs for the finite difference numerical derivatives
'''
import numpy as np
import math

def coef_r1 (n):
    A=np.zeros((5,5),dtype='float')
    C=np.zeros(5, dtype='float')
    C[n]=1

    for i in range(5):
        A[:,i]=i
        for j in range(5):
            A[j,i]=(i**j)/math.factorial(j)
    
    return np.matmul( np.linalg.inv(A) ,C)


def coef_r1_npoint (n, nd):
    A=np.zeros((nd,nd),dtype='float')
    C=np.zeros(nd, dtype='float')
    C[n]=1

    for i in range(nd):
        A[:,i]=i
        for j in range(nd):
            A[j,i]=(i**j)/math.factorial(j)
    
    return np.matmul( np.linalg.inv(A) ,C)

print(coef_r1_npoint (2, 5))


print(coef_r1_npoint (2, 7))

def coef_symmetric_center (n, nd):
    A=np.zeros((nd,nd),dtype='float')
    C=np.zeros(nd, dtype='float')
    C[n]=1

    for i in range(nd):
        A[:,i]=i
        for j in range(nd):
            if (j>0):
                A[j,i]=(i-2)**j
            A[j,i]=A[j,i]/math.factorial(j)
    
    return np.matmul( np.linalg.inv(A) ,C)
    
print( coef_symmetric_center (1,5)  )    
