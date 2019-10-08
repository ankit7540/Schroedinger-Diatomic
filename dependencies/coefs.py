'''
Programs for the coefs for the finite difference numerical derivatives
'''
import numpy as np
import math
import timeit
import functools
#-------------------------------------------------------------
def coef_symmetric_center (n, nd):
    '''
    finite difference coef for symmetric (center) derivative
    Params:
        n   = nth order derivative, non-negative integer
        nd  = nd point stencil, non-negative integer
        condition : n < nd 
    
    Returns:
        array[float64] of length =nd
    '''

    A=np.zeros((nd,nd))
    C=np.zeros(nd)
    C[n]=1

    for i in range(nd):
        A[:,i]=i
        for j in range(nd):
            if (j>0):
                A[j,i]=(i-(math.floor(nd/2)))**j
            else:
                A[j,i]=i**j    
            A[j,i]=A[j,i]/math.factorial(j)

    return np.linalg.solve(A, C)
#-------------------------------------------------------------

def coef_asymmetric_forward (n, nd):
    '''
    finite difference coef for asymmetric (forward) derivative
    Params:
        n   = nth order derivative, non-negative integer
        nd  = nd point stencil, non-negative integer
        condition : n < nd 
    
    Returns:
        array[float64] of length =nd
    '''

    A=np.zeros((5,5),dtype='float')
    C=np.zeros(5, dtype='float')
    C[n]=1

    for i in range(5):
        A[:,i]=i
        for j in range(5):
            A[j,i]=(i**j)/math.factorial(j)
    
    return np.matmul( np.linalg.inv(A) ,C)

#-------------------------------------------------------------

def coef_r1_npoint (n, nd):
    A=np.zeros((nd,nd),dtype='float')
    C=np.zeros(nd, dtype='float')
    C[n]=1

    for i in range(nd):
        A[:,i]=i
        for j in range(nd):
            A[j,i]=(i**j)/math.factorial(j)
    
    return np.linalg.solve(A, C)

#-------------------------------------------------------------


#----------------------------------
# timing the functions 
A=1
B=13
t = timeit.Timer(functools.partial(coef_symmetric_center, A, B)) 
print (t.timeit(10))

#----------------------------------


    

#X=( coef_symmetric_center (1,11)  )    
#print(X)