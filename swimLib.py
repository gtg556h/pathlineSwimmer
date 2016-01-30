import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg
import matplotlib.animation as animation

def genWaveform1(t,omega):
    # Standard waveform used in NatureComm paper

    m = np.sin(omega/2.0*t)
    m = (np.abs(m) - 0.3)/0.7
    m[np.where(m<0)[0]] = 0

    return m

def genWaveform2(t,omega):
    # Pefect sinusoidal, 2-sided driving function

    m = np.sin(omega*t)
    #m = (np.abs(m) - 0.3)/0.7
    #m[np.where(m<0)[0]] = 0

    return m

def genWaveform3(t,omega):
    # Nominally tweaked waveform form conference proceedings etc...

    m = np.sin(omega/2.0*t)
    m = (np.abs(m) - 0.25)/0.75
    m[np.where(m<0)[0]] = 0

    return m



def calc_x_s(x, ds):
    return 0

def calc_x_ss(x, ds):
    return 0

def calc_x_sss(x, ds):
    return 0

def calc_x_ssss(x, ds):
    return 0

def calcNormals(x, ds):
    # e_t = (dx/ds, dy/ds)
    # e_n = (-dy/ds, dx/ds)
    x_s = calc_x_s(x,s)
    return 0











