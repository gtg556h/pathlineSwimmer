import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg
import matplotlib.animation as animation
import pdb
import numDiffLib as ndl


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

def gen_x_s(ns, ds):
    A = ndl.genOperator(ns, 1, ds)
    return lambda x: np.dot(A,x.T).T

def gen_x_ss(ns, ds):
    A = ndl.genOperator(ns, 2, ds)
    return lambda x: np.dot(A, x.T).T

def gen_x_sss(ns, ds):
    A = ndl.genOperator(ns, 3, ds)
    return lambda x: np.dot(A, x.T).T

def gen_x_ssss(ns, ds):
    A = ndl.genOperator(ns, 4, ds)
    return lambda x: np.dot(A, x.T).T

#############################











