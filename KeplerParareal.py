import sys 
sys.path.append('Euler')
import parareal as para
import RK2
import RK4
import HamiltonianFunctions as func
import numpy as np
import matplotlib.pyplot as plt
##Constants
T = 50
N = 50
delta_T = T/N
e = 0.7
a =1
b = np.sqrt(1-e**2)
d = 1-e**2
y0 = np.array([0,np.sqrt((1+ e)/(1-e)),1-e,0])
dtg = delta_T/2
dtf = 0.1
epsilon = 0.001

## Computations
tab_y, tab_err = para.parareal(func.Kepler, RK4.RK4, RK2.RK2,T, y0, dtf,dtg, N, epsilon)
plt.plot(tab_y[2], tab_y[3],'o',label='parareal')
tab_t = np.linspace(0,T,N+1)
sol_ex = np.array([a*np.cos(tab_t)-e,b*np.sin(tab_t)])
plt.plot(sol_ex[0],sol_ex[1],'x',label='exact')
plt.legend()
plt.show()

plt.plot(np.reshape(np.array([range(N)]),N),tab_err)
plt.yscale('log')
plt.show()
