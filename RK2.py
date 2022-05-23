##Imports
import numpy as np
import matplotlib.pyplot as plt

#Functions
def y_n(yn, tn,f,dt):
    return yn + dt * f(tn, yn)/2

def dy_n(yn,tn,f,dt):
    return f(tn + dt/2, y_n(yn, tn, f, dt))

def small_RK2(f,t0,y0,dt):
    return y0 + dt * dy_n(y0, t0 + dt, f, dt)

def tab_RK2(f,t0,T,y0,dt):
    y = y0
    m = int((T-t0)/dt)
    tab_y = np.zeros((len(y0), m+1))
    for n in range(m+1):
        tab_y[:,n] = y
        y = small_RK2(f, t0 + (n-1) * dt, y, dt)
    return tab_y

def RK2(f,t0,T,y0,dt):
    y = y0
    m = int((T-t0)/dt)
    for n in range(m+1):
        y = small_RK2(f, t0 + (n-1) * dt, y, dt)
    return y