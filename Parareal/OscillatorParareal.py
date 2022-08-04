import sys
sys.path.append('Euler')
import EulerImplicitOscillator as EIO
import parareal as para
import RK2
import RK4
import StormerOscillator as SO
import HamiltonianFunctions as func
import numpy as np
import matplotlib.pyplot as plt



# Constants
T = 6
N = 60  
y0 = np.array([1, 0])
delta_t = T/N
dtg = delta_t
dtf = dtg/100
epsilon = 0.001

# Computations
# tab_y = para.parareal(func.Oscillator, RK4.RK4,
#                       SO.Stormer, T, y0, dtf, dtg, N, epsilon)

# plt.plot(tab_y[0, :], tab_y[1, :], 'o', label='parareal')
# tab_t = np.linspace(0, T, N+1)

# sol_ex = np.array([np.cos(tab_t), np.sin(tab_t)])
# plt.plot(sol_ex[0], sol_ex[1], 'x', label='exact')
# plt.legend()
# plt.show()

##Error computation
#Error on the solution
tab_err = para.parareal_err(func.Oscillator, RK4.tab_RK4, RK4.tab_RK4, T, y0, dtf, dtg, N, epsilon)
plt.loglog(np.linspace(0,T,N+1),tab_err)
plt.show()

#Error on the Hamiltonian

# tab_Ha = para.parareal_Ha_Oscillator(func.Oscillator, RK4.RK4, RK2.RK2,T, y0, dtf,dtg, N, epsilon)

# for k in range(1,6):
# 	plt.plot(np.linspace(0,T,N+1)[:], (np.abs(tab_Ha[k]-0.5*np.ones(len(tab_Ha[k]))))[:],label=f"k={k}")
# 	plt.yscale('log')

# print("dtf", dtf)
# sol_RK4= RK4.RK4(func.Oscillator,0,T,y0,dtf)
# sol_RK2 = RK2.RK2(func.Oscillator,0,T,y0,dtg)
# Ha_RK4 = []
# Ha_RK2 = []
# print("RK4",np.shape(sol_RK4))
# print("RK2",np.shape(sol_RK2))
# for k in range(np.shape(sol_RK2)[1]):
# 	Ha_RK4.append(0.5*(sol_RK4[0,100*k]**2 +sol_RK4[1,100*k]**2))

# for k in range(np.shape(sol_RK2)[1]):
# 	Ha_RK2.append(0.5*(sol_RK2[0,k]**2 +sol_RK2[1,k]**2))
# print("HA4",np.shape(Ha_RK4))
# print("HA2",np.shape(Ha_RK2))
# plt.plot(np.linspace(0,T,int(T/dtg)+1), (Ha_RK4-Ha_RK4[0]*np.ones(len(Ha_RK4))),label='RK4')
# plt.plot(np.linspace(0,T,int(T/dtg)+1),(Ha_RK2-Ha_RK2[0]*np.ones(len(Ha_RK2))),label='RK2')
# plt.legend()
# plt.show()
# tab_err = []
# delta_t = 0.1
# dtf =0.1
# dtg=0.1
# for m in range(10,101):
# 	tab_y = para.parareal(func.Oscillator, RK4.RK4, SO.Stormer,delta_t*m, y0, dtf,dtg, m, epsilon)
# 	sol_ex = np.array([np.cos(np.linspace(0,delta_t*m,m+1)),np.sin(np.linspace(0,delta_t*m,m+1))])
# 	tab_err.append(np.linalg.norm(tab_y[:,10]-sol_ex[:,10]))

# plt.loglog(np.linspace(10,100,91),tab_err,label='error')
# plt.legend()
# plt.title('Error of approximate')
# plt.show()
