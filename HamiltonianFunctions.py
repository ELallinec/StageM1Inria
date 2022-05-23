import numpy as np

def Oscillator(tn,yn):
    return np.array([-yn[1],yn[0]])

def Kepler(t,y):
    p = y[0:2]
    q = y[2:4]
    const = (q[0]**2 +q[1]**2)**(-3/2)
    return np.reshape(np.array([-q*const,p]),4)
