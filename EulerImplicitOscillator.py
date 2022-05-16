#Euler implicit symplectic

import numpy as np 
import matplotlib.pyplot as plt

h = 0.1
m = 100
T = h*m
p0 = 0
q0 = 0
#Euler for H(p,q) = 1/2 * (p² +q²), grad_p(H)= p, grad_q(H)=q
def Euler_harmonic(p,q,m,h):
    for n in range(m):
        p = p -h*q
        q = q -h*p
        
    return np.array([p,q])

[p,q] = Euler_harmonic(p0,q0,m,h) 

plt.plot(p,q)