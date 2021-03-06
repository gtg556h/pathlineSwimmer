import numpy as np
import scipy as sci
import scipy.linalg
import math
import pdb




def returnDerivative(nTerms, position, order):
    # position is zero indexed!

    if (position+1) > nTerms:
        print("yo, position is outside of nTerms!!!")

    A = np.zeros((nTerms, nTerms))

    for ii in range(0,A.shape[0]):
        for jj in range(0,A.shape[0]):
            A[ii,jj] = (jj-position)**(ii)/float(math.factorial(math.fabs(ii)))
    b = np.zeros((nTerms,1))
    b[order] = 1

    c = scipy.linalg.solve(A,b)
    return c.T
    



def genOperator(n, order, ds=1):
    # n is size of square operator
    # order is order of differentiation

    A = np.zeros([n, n])

    specialCases = np.int(np.floor((order+1)/2))

    if np.mod(order, 2) == 0:
        centralTerms = order + 1
    else:
        centralTerms = order + 2

    for j in range(0, specialCases):
        #pdb.set_trace()
        A[j,0:order+2] = returnDerivative(order+2, j, order)

    for j in range(specialCases, n - specialCases):
        #pdb.set_trace()
        A[j, j-specialCases:j+specialCases+1] = returnDerivative(centralTerms, np.int(np.floor(centralTerms/2)), order)

    for j in range(n-specialCases,n):
        #pdb.set_trace()
        A[j, -order-2:] = returnDerivative(order+2, j-n+2+order, order)
        
    return A / ds**order
