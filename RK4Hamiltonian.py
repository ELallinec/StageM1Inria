# -*- coding: utf-8 -*-
#%%
"""
Created on Thu Mar 31 17:05:36 2022

@author: ewen7
"""

import numpy as np
import matplotlib.pyplot as plt


##
h = 0.01
m = 100

def k1(f,tn,yn):
    return f(tn,yn)

def k2(f,tn,yn,h):
    return f(tn+h/2,yn+h*k1(f,tn,yn)/2)

def k3(f,tn,yn,h):
    return f(tn+h/2,yn+h*k2(f,tn,yn,h)/2)

def k4(f,tn,yn,h):
    return f(tn+h,yn+h*k3(f,tn,yn,h))

##f1 = J^-1 Grad(H(p,q)) avec H(p,q) = 1/2(p^2+q^2)
def f1(tn,yn):
    try:
        n = len(yn)//2
        return np.hstack([[-yn[0]],[yn[1]]])
    except TypeError:
        print("Le type de donnee n'est pas bon ")
        pass
    

def RK4(f,t0,y0,h,m):
    y = y0
    t = t0
    tab_y = np.array([y])
    tab_t = np.array([t])
    for k in range(m):
        y = y + 1/6 * h * ( k1(f,t,y) + 2 * k2(f,t,y,h) + 2 * k3(f,t,y,h) + k4(f,t,y,h) )   
        t += h
        tab_y = np.concatenate((tab_y,[y]))
        tab_t= np.concatenate((tab_t,[t]))
    return tab_y, tab_t

[tab_y,tab_t] = RK4(f1,0,np.array([0,0]),h,m)

#print(tab_y)
plt.plot(tab_y[:,0],tab_y[:,1])

# %%
