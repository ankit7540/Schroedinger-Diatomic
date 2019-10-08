'''
Programs for the coefs for the finite difference numerical derivatives
'''

def coef_r1 (n):
    A=np.zeros((5,5))
    C=np.zeros(5)
    
    for i in range(5):
        A[:,i]=i
        for j in range(5):
            A[j,i]=(i^j)/factorial(j)
    
    print(A)        