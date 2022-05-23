#Euler implicit symplectic

#Imports
import numpy as np 
import matplotlib.pyplot as plt

#Constants
p0 = 1
q0 = 0
dt = 0.1
n = 1000
T = 1

#Euler for H(p,q) = 1/2 * (p² +q²), grad_p(H)= p, grad_q(H)=q
def Euler_Symplectic_harmonic(p,q,m,h):
    tab = np.zeros((2,m))
    Ha = np.zeros(m)
    for n in range(m):
        Ha[n] = (p**2 + q**2)/2
        tab[0,n] = p
        tab[1,n] = q
        p = p -h*q
        q = q +h*p
    return tab,Ha

#Plot approx vs exact solution
tab,Ha= Euler_Symplectic_harmonic(p0,q0,n,dt)
plt.plot(tab[0],tab[1])
plt.title("Euler Implicit")
plt.xlabel("Position") 
plt.ylabel("Velocity")
plt.show()


x_ex = np.cos(np.linspace(0,n*dt,1000))
v_ex = np.sin(np.linspace(0,n*dt,1000))
plt.plot(x_ex,v_ex)
plt.title("Exact solution")
plt.xlabel("Position") 
plt.ylabel("Velocity")
plt.show()

#Hamiltonian through time
plt.plot(np.linspace(0,n*dt,1000),Ha)
plt.xlabel("Time")
plt.ylabel("Hamiltonian")
plt.title("Hamiltonian evolution through time")
plt.show()


#Error computation
tab_err = np.zeros((2,10))
for p in range(6, 16):
    h = T/(2**p)
    tab,Ha= Euler_Symplectic_harmonic(p0,q0,(2**p),h) 
    x_ex = np.cos(np.linspace(0,T,(2**p)))
    v_ex = np.sin(np.linspace(0,T,(2**p)))
    err_x = np.max(np.abs(np.subtract(tab[0],x_ex)))
    err_v = np.max(np.abs(np.subtract(tab[1],v_ex)))
    tab_err[0,p-6] = err_x
    tab_err[1,p-6] = err_v

#Order computation
slope_x, intercept_x = np.polyfit(np.log([2**p for p in range(6,16)]),np.log(tab_err[0]),1)
slope_v, intercept_v = np.polyfit(np.log([2**p for p in range(6,16)]),np.log(tab_err[1]),1)

#Error plot  
plt.loglog([2**p for p in range(6,16)],tab_err[0],label=f"pos_err, order: {-slope_x}")
plt.loglog([2**p for p in range(6,16)],tab_err[1],label=f"v_err, order: {-slope_v}")
plt.title("Errors in function of m")
plt.legend()
plt.show()

