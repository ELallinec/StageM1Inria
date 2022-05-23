#Imports
import numpy as np 
import matplotlib.pyplot as plt

#Constants
p0 = 1
q0 = 0
dt = 0.1
n = 1000
T = 1

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


def Euler_harmonic_imp(p,q,m,h):
    tab = np.zeros((2,m))
    Ha = np.zeros(m)
    for n in range(m):
        Ha[n] = (p**2 + q**2)/2
        tab[0,n] = p
        tab[1,n] = q
        q = (q+h*p)/(1+h**2)
        p = p -h*q

    return tab, Ha

def Euler_harmonic_ex(p,q,m,h):
    tab = np.zeros((2,m))
    Ha = np.zeros(m)
    for n in range(m):
        Ha[n] = (p**2 + q**2)/2
        tab[0,n] = p
        tab[1,n] = q
        temp =p
        p = p -h*q
        
        q = q +h*temp
    return tab,Ha

tab_imp, Ha_imp = Euler_harmonic_imp(p0,q0,n,dt)
tab_ex, Ha_ex = Euler_harmonic_ex(p0,q0,n,dt)
tab_symp, Ha_symp = Euler_Symplectic_harmonic(p0,q0,n,dt)


tab_t = np.linspace(0,n*dt,n)
plt.plot(tab_t, Ha_imp/np.max(Ha_imp), label='Implicit')
plt.plot(tab_t, Ha_ex/np.max(Ha_ex), label ='Explicit')
plt.plot(tab_t, Ha_symp/np.max(Ha_symp), label = 'Symplectic')
plt.legend()
plt.title("Comparison of Hamiltonian through time using \n Euler implicit, explicit and implicit symplectic")
plt.show()