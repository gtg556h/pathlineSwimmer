import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg
import matplotlib.animation as animation
import pdb



#####################################################
#####################################################


def genWaveform1(t,omega):
    # Standard waveform used in NatureComm paper

    m = np.sin(omega/2.0*t)
    m = (np.abs(m) - 0.3)/0.7
    m[np.where(m<0)[0]] = 0

    return m


#############################

def genWaveform2(t,omega):
    # Pefect sinusoidal, 2-sided driving function

    m = np.sin(omega*t)
    #m = (np.abs(m) - 0.3)/0.7
    #m[np.where(m<0)[0]] = 0

    return m

#############################

def genWaveform3(t,omega):
    # Nominally tweaked waveform form conference proceedings etc...

    m = np.sin(omega/2.0*t)
    m = (np.abs(m) - 0.25)/0.75
    m[np.where(m<0)[0]] = 0

    return m


#############################


####################################################
####################################################
####################################################

# REWRITE THE FOLLOWING FUNCTIONS TO GENERATE LAMBDA FUNCTIONS
# THAT THEN GET MAPPED TO CLASS METHODS

def calc_x_s(x, ds):
    # Argument:  x vector at *single time point*
    # with dimensions [2, ns]
    ns = x.shape[1]

    A = np.zeros([ns,ns])
    A[0, 0:2] = np.array([-1, 1])/ds
    A[-1,-2:] = np.array([-1,1])/ds
    for j in range(1,A.shape[0] - 1):
        A[j,j-1:j+2] = np.array([-1,0,1])/(2*ds)
    
    x_s = np.dot(A, x.T).T
    return x_s


#############################

def calc_x_ss(x, ds):

    ns = x.shape[1]

    A = np.zeros([ns,ns])
    A[0,0:4] = np.array([2, -5, 4, -1])/ds**2
    A[-1,-4:] = np.array([-1, 4, -5, 2])/ds**2
    for j in range(1, A.shape[0] - 1):
        A[j,j-1:j+2] = np.array([1, -2, 1])/ds**2

    x_ss = np.dot(A, x.T).T

    return x_ss
    

#############################

def calc_x_sss(x, ds):

    ns = x.shape[1]

    A = np.zeros([ns,ns])
    A[0,0:5] = np.array([-2.5, 9, -12, 7, -1.5])/ds**3
    A[1,0:5] = np.array([-1.5, 5, -6, 3, -0.5])/ds**3
    A[-2,-5:] = np.array([0.5, -3, 6, -5, 1.5])/ds**3
    A[-1,-5:] = np.array([1.5, -7, 12, -9, 2.5])/ds**3
    for j in range(2, A.shape[0] - 2):
        A[j,j-2:j+3] = np.array([-0.5, 1, 0, -1, 0.5])/ds**3

    x_sss = np.dot(A, x.T).T

    return x_sss


#############################

def calc_x_ssss(x, ds):

    ns = x.shape[1]
    A = np.zeros([ns, ns])
    A[0,0:6] = np.array([3, -14, 26, -24, 11, -2])/ds**4
    A[1,0:6] = np.array([2, -9, 16, -14, 6, -1])/ds**4
    A[-2,-6:] = np.array([-1, 6, -14, 16, -9, 2])/ds**4
    A[-1,-6:] = np.array([-2, 11, -24, 26, -14, 3])/ds**4
    for j in range(2, A.shape[0] - 2):
        A[j,j-2:j+3] = np.array([-0.5, 1, 0, -1, 0.5])/ds**4

    x_ssss = np.dot(A, x.T).T

    return x_ssss


#############################

def calcNormals(x, ds):
    # e_t = (dx/ds, dy/ds)
    # e_n = (-dy/ds, dx/ds)
    x_s = calc_x_s(x,ds)

    e_t = np.zeros_like(x)
    e_n = np.zeros_like(x)

    e_t[0,:] = x_s[0,:]
    e_t[1,:] = x_s[1,:]

    e_n[0,:] = -x_s[1,:]
    e_n[1,:] = x_s[0,:]

    return e_t, e_n
    

###############################










