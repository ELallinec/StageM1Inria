##Imports
import numpy as np
import matplotlib.pyplot as plt

##Functions
def tab_Stormer(f,t0,T,y0,dt):
    m = int((T-t0)/dt)
    tab = np.zeros((len(y0), m+1))
    p = y0[:2]
    q = y0[2:]
    for k in range(m+1):
        tab[:,k] = np.array([p[0],p[1], q[0],q[1]])
        p_temp = p - (dt / 2) *q * ((q[0]**2+q[1]**2)**(-1.5))
        q = q + dt * p_temp
        p = p_temp -(dt / 2) * q * ((q[0]**2+q[1]**2)**(-1.5))
    return tab

def Stormer(f,t0,T,y0,dt):
    m = int((T-t0)/dt)
    p = y0[0:2]
    q = y0[2:]
    for k in range(m+1):
        p_temp = p - (dt / 2) * q * ((q[0]**2+q[1]**2)**(-1.5))
        q = q + dt * p_temp
        p = p_temp -(dt / 2) * q *((q[0]**2+q[1]**2)**(-1.5))
    return np.hstack((p,q))