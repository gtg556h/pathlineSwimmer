import numpy as np
import matplotlib.pyplot as plt
import swimLib as sl
import flagellaClass as fc 

# x[i,s,t], i=[1,2] for x, y coord

L = 2000
ds = 3.0

s = np.arange(0, L+ds, ds)
x0 = np.zeros([2,s.shape[0]])
x0[0,:] = s.copy()
x0[1,:] = np.zeros_like(s)


#omega = 3.6*2*np.pi
omega = 3.25*2*np.pi
nT = 4.0  # n Periods
dt = 0.001

drivingFunction = 3

E = 3.86
I = 7.9**3*22.0/12.0
A = E*I



zetaN = 9.25E-9
zetaT = zetaN / 2.0

moment = 37.37
mStart = 445.0   # unit length, not grid number!
mEnd = mStart + 60.0

BC = [0,0]
twist = [0,0]
shift = [0,0]
shiftAmp = 1
twistAmp = 1

params = {'L':L, 'ds':ds, 's':s, 'x0':x0, 'omega':omega, 'nT':nT, 'dt':dt, 'E':E, 'A':A,  'zetaN':zetaN, 'zetaT':zetaT, 'moment':moment, 'mStart':mStart, 'mEnd':mEnd, 'BC':BC, 'twist':twist, 'shift':shift, 'shiftAmp':shiftAmp, 'twistAmp':twistAmp, 'drivingFunction':drivingFunction}  


s = fc.flagella(params)
#s1.actuator()
#s1.numSolve()
#s1.propulsionCalc()
#print(np.mean(s1.Ux))
#s1.plotDisp(16,1)




