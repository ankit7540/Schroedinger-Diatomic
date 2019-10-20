'''
Programs for the coefs for the finite difference numerical derivatives
'''
import numpy as np
import math

#import timeit
#import functools
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

def coef_forward (n, nd ):
    '''
    finite difference coef for forward derivative
    Params:
        n   = nth order derivative, non-negative integer
        nd  = nd point stencil, non-negative integer
        condition : n < nd 
    
    Returns:
        array[float64] of length =nd
    '''

    A=np.zeros((nd,nd),dtype='float')
    C=np.zeros(nd, dtype='float')
    C[n]=1

    for i in range(nd):
        A[:,i]=i
        for j in range(nd):  
            A[j,i]=(i**j)/math.factorial(j)
    #print(A)        
            
    return np.linalg.solve(A, C)
    
#-------------------------------------------------------------

def coef_backward (n, nd ):
    '''
    finite difference coef for backward derivative
    Params:
        n   = nth order derivative, non-negative integer
        nd  = nd point stencil, non-negative integer
        condition : n < nd 
    
    Returns:
        array[float64] of length =nd
    '''

    A=np.zeros((nd,nd),dtype='float')
    C=np.zeros(nd, dtype='float')
    C[n]=1

    for i in range(nd):
        A[:,i]=i
        for j in range(nd):  
            A[j,i]=(i**j)/math.factorial(j)
        
            if j % 2 != 0:
                A[j,i]=-1*A[j,i] 
        
    #print(A)            
    return np.linalg.solve(A, C)
    
#-------------------------------------------------------------
    

def coef_forward_asymmetric (n, nd, m ):
    '''
    finite difference coef for symmetric (center) derivative
    Params:
        n   = nth order derivative, non-negative integer
        nd  = nd point stencil, non-negative integer
        condition : n < nd 
        m   = 1 : 1 point backward, other forward
              2 : 2 point backward, other forward
              .
              M : M point backward, other forward
    
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
                A[j,i]=(i-m)**j
            else:
                A[j,i]=i**j    
            A[j,i]=A[j,i]/math.factorial(j)

    return np.linalg.solve(A, C)
    
#-------------------------------------------------------------

def coef_backward_asymmetric (n, nd, m ):
    '''
    finite difference coef for symmetric (center) derivative
    Params:
        n   = nth order derivative, non-negative integer
        nd  = nd point stencil, non-negative integer
        condition : n < nd
        m   = 1 : 1 point backward, other forward
              2 : 2 point backward, other forward
              .
              M : M point backward, other forward        
    
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
                A[j,i]=(i-m)**j
            else:
                A[j,i]=i**j    
            A[j,i]=A[j,i]/math.factorial(j)
            
            if j % 2 != 0:
                A[j,i]=-1*A[j,i] 

    return np.linalg.solve(A, C)    

#-------------------------------------------------------------


#----------------------------------
# timing the functions 
#A=1
#B=13
#M=1
#t = timeit.Timer(functools.partial(coef_asymmetric_forward, A, B, M)) 
#print (t.timeit(10))

#----------------------------------


#print(coef_symmetric_center( A, B) )
S1=coef_backward( 1, 5 ) 

    

#X=( coef_symmetric_center (1,11)  )    
#print(X)
S2=coef_backward_asymmetric (1, 5, 1 ) 
S3=coef_backward (1, 5 )
print(S2)
print(S3)

print ("\n")

S2=coef_backward_asymmetric (2, 5, 1 ) 
S3=coef_backward (2, 5 )
print(S2)
print(S3)

