#Imports
import numpy as np 
import matplotlib.pyplot as plt

#Approx
def tab_Euler(f,t0,T,y0,h):
    m = int((T - t0) / h)
    tab = np.zeros((2, m+1))   
    p = y0[0]
    q = y0[1]
    for n in range(m+1):
        tab[0,n] = p
        tab[1,n] = q
        q = (q + h *  p) / (1 + h**2)
        p = p -h*q
    return tab

def Euler(f,t0,T,y0,h):
    m = int((T - t0) / h)
    p = y0[0]
    q = y0[1]
    for n in range(m+1):
        q = (q + h * p) / (1 + h**2)
        p = p - h * q
    return np.array([p,q])