
##Imports
import numpy as np
import matplotlib.pyplot as plt

##Functions
def tab_Stormer(f,t0,T,y0,dt):
    m = int((T-t0)/dt)
    tab = np.zeros((2, m+1))
    p = y0[0]
    q = y0[1]
    for k in range(m+1):
        tab[:,k] = np.array([p, q])
        p_temp = p - (dt / 2) *q
        q = q + dt * p_temp
        p = p_temp -(dt / 2) * q
    return tab

def Stormer(f,t0,T,y0,dt):
    m = int((T-t0)/dt)
    p = y0[0]
    q = y0[1]
    for k in range(m+1):
        p_temp = p - (dt / 2) * q
        q = q + dt * p_temp
        p = p_temp -(dt / 2) * q
    return np.array([p, q])