##Imports
import numpy as np
import matplotlib.pyplot as plt


##Constants
dt = 0.001
n = 10000
p0 = 1
q0 = 0

##Functions
def p_nhalf(p,q,h):
    return p -h*q/2

def q_n(p,q,h):
    return q+ h/2*p_nhalf(p,q,h)

def p_n(p,q,h): 
    return p_nhalf(p,q,h) -h*q_n(p,q,h)/2

def Stormer(p0,q0,m,h):
    tab = np.zeros((2,m))
    Ha = np.zeros(m)
    p = p0
    q = q0
    for k in range(m):
        Ha[k] = (p**2 + q**2)/2
        tab[0,k] = p
        tab[1,k] = q
        
        q = q_n(p,q,h)
        p = p_n(p,q,h)
    return tab,Ha

tab, Ha= Stormer(p0,q0,n,dt)

#First plot
plt.plot(tab[0],tab[1])
plt.title("Stormer-Verlet")
plt.xlabel("Position") 
plt.ylabel("Velocity")
plt.show()

v_ex = np.sin(np.linspace(0,n*dt,n))
x_ex = np.cos(np.linspace(0,n*dt,n))
plt.plot(x_ex,v_ex)
plt.title("Exact solution")
plt.xlabel("Position") 
plt.ylabel("Velocity")
plt.show()

#Hamiltonian through time
plt.plot(np.linspace(0,n*dt,n),Ha)
plt.xlabel("Time")
plt.ylabel("Hamiltonian")
plt.title("Hamiltonian evolution through time")
plt.show()
