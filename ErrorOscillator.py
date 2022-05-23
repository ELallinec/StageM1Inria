#!/usr/bin/env python
# coding: utf-8

# # Harmonic oscillator 

# We have an Hamiltonian system
# $$\dot{y} = J^{-1} \nabla H(y) $$
# with $y = (p,q) \in \mathbb{R}^{2d}$ and $J=\begin{pmatrix}0 & I \\ -I & 0 \end{pmatrix} \in M_{2d}(\mathbb{R})$

# In this case  $d=1$ and we have the following hamiltonian
# \begin{equation}
# H(p,q) = \frac{1}{2}\left(p^2 + q ^2\right)
# \end{equation}

# Thus
# 
# $$
# \left\{\begin{matrix}
#        \frac{\partial p}{\partial t} &=& - q \\
#        \frac{\partial q}{\partial t} &=& p\\
# \end{matrix}\right.
# $$
# 

# ### Imports

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
#matplotlib.rcParams['figure.figsize'] = [10, 5]
def parareal(F, f, g, T, y0, dtf, dtg, N, epsilon):
    delta_t = T/N
    kmax = N
    lambda_tab_old = np.zeros((len(y0), N+1))
    lambda_tab = np.zeros((len(y0), N+1))
    lambda_f = np.zeros((len(y0),N+1))
    
    
    lambda_tab_old[:, 0] = y0
    for n in range(1, N+1):
        lambda_tab_old[:, n] = g(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtg)
    
    k = 1
    lambda_f[:,0] = y0
    lambda_tab[:, 0] = y0
    while k <= kmax:
        
        
        for n in range(1,N+1):
            lambda_f[:, n] = f(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtf) 
        
        
        for n in range(1, N+1):
            lambda_tab[:, n] = lambda_f[:,n]- g(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtg)
            lambda_tab[:, n] = lambda_tab[:, n] + g(F, (n-1)*delta_t, n*delta_t, lambda_tab[:, n-1], dtg)
        
        lambda_tab_old = np.copy(lambda_tab)
        
        print("k:", k)
        k += 1
    return lambda_tab


def parareal_err(F, f, g, T, y0, dtf, dtg, N, epsilon):
    delta_t = T/N
    kmax = N
    tab_err = []
    lambda_tab_old = np.zeros((len(y0), N+1))
    lambda_tab = np.zeros((len(y0), N+1))
    sol_ex = f(F, 0, T, y0, delta_t)
    lambda_tab_old[:, 0] = y0
    for n in range(1, N+1):
        lambda_tab_old[:, n] = g(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtg)[:,-1]
    k = 1
    tab_err.append(np.linalg.norm(sol_ex-lambda_tab_old))

    while k <= kmax:
        for n in range(1, N+1):
            lambda_tab[:, 0] = y0
            lambda_tab[:, n] = f(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtf)[:,-1] - g(
                F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtg)[:,-1]
            lambda_tab[:, n] = lambda_tab[:, n] + \
                g(F, (n-1)*delta_t, n*delta_t, lambda_tab[:, n-1], dtg)[:,-1]
        lambda_tab_old = np.copy(lambda_tab)
        print("k:", k)
        
        tab_err.append(np.linalg.norm(sol_ex-lambda_tab))
        if np.abs(tab_err[k]-tab_err[k-1])<epsilon:
            break
        k += 1
    return np.array(tab_err)



def parareal_Ha_Oscillator(F, f, g, T, y0, dtf, dtg, N, epsilon):
    delta_t = T/N
    kmax = 6
    lambda_tab_old = np.zeros((len(y0), N+1))
    lambda_f = np.zeros((len(y0),N+1))
    lambda_tab = np.zeros((len(y0), N+1))
    tab_Ha = 0.5 * np.ones((kmax+2, N+1))
    tab_f = np.zeros((len(y0),N+1))

    tab_f[:,0] = y0
    tab_Ha[0,0]= 0
    for n in range(1,N+1):
        tab_f[:, n] = f(F, (n-1)*delta_t, n*delta_t, tab_f[:, n-1], dtf)
        tab_Ha[0, n] = tab_Ha[0, n] - 0.5 * (tab_f[0, n]**2 + tab_f[1, n]**2)
    
    lambda_tab_old[:, 0] = y0
    tab_Ha[1,0] = 0
    for n in range(1, N+1):
        lambda_tab_old[:, n] = g(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtg)
        tab_Ha[1, n] = tab_Ha[1, n] - 0.5 * (lambda_tab_old[0, n]**2 + lambda_tab_old[1, n]**2)
    
    k = 1
    while k <= kmax:
        tab_Ha[k+1,0] = 0
        lambda_f[:,0] = y0

        for n in range(1,N+1):
            lambda_f[:, n] = f(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtf) 
        
        lambda_tab[:, 0] = y0
        for n in range(1, N+1):
            lambda_tab[:,n] = lambda_f[:,n] - g(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtg)
            lambda_tab[:, n] = lambda_tab[:, n] + g(F, (n-1)*delta_t, n*delta_t, lambda_tab[:, n-1], dtg)
            tab_Ha[k+1, n] = tab_Ha[k+1, n] -0.5 * (lambda_tab[0, n]**2 + lambda_tab[1, n]**2)
        lambda_tab_old = np.copy(lambda_tab)
        
        print("k:", k)
        k += 1
    return tab_Ha



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
        p_temp = p - (dt/2)*q
        q = q + dt * p_temp
        p = p_temp -(dt / 2) * q
    return np.array([p, q])

def Oscillator(tn,yn):
    return np.array([-yn[1],yn[0]])
# ### Constants

# In[2]:


T = 1000
N = 10000
y0 = np.array([0, 1])
delta_t = T/N
dtg = delta_t
dtf = dtg/100
epsilon = 1/N


# ## Plot

tab_Ha = parareal_Ha_Oscillator(Oscillator, Stormer, Stormer, T, y0, dtf, dtg, N, epsilon)

kmax = len(tab_Ha)

tab_t = np.linspace(0,T,N+1)

tab_Ha0 = 0.5 * np.ones(len(tab_Ha[0,:]))

plt.loglog(tab_t[1:],np.abs(tab_Ha[0,1:]), label ="f solver")
plt.loglog(tab_t[1:],np.abs(tab_Ha[1,1:]), label ="g solver")

for k in range(2,kmax):
    plt.loglog(tab_t[1:],np.abs(tab_Ha[k,1:]), label =f"k={k-1}")
    plt.xlabel('Time')
    plt.ylabel('Error')

plt.legend()
plt.savefig('Error.png')


# In[7]:


for k in range(kmax):
    slope, intercept = np.polyfit(np.log(tab_t)[1:],np.log(np.abs(tab_Ha[k,1:])),1)
    print(f"pente {k-1}", slope)



# In[ ]:




