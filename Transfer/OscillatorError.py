import StormerOscillator as SO
import HamiltonianFunctions as func
import numpy as np
import matplotlib.pyplot as plt
import parareal as para

T = 100000
N = 1000000
y0 = np.array([0, 1])
delta_t = T/N
dtg = delta_t
dtf = delta_t/100
epsilon = 1/N
 

tab_y = para.parareal_bis(func.Oscillator, SO.Stormer,SO.Stormer, T, y0, dtf, dtg, N-1, epsilon)

tab_t = np.linspace(0, T, N+1)
sol_ex = np.array([-np.sin(tab_t),np.cos(tab_t)])


tab_Ha0 = 0.5 * np.ones(len(tab_y[0]))
tab_Ha = para.Ha_err(tab_y,tab_Ha0,func.HaOscillator)
kmax = len(tab_Ha)


plt.loglog(tab_t,tab_Ha[0], label ="f solver")
plt.loglog(tab_t,tab_Ha[1], label ="g solver")

for k in range(2,kmax):
    plt.loglog(tab_t,tab_Ha[k], label =f"k={k-1}")
    plt.xlabel('Time')
    plt.ylabel('Error')
    plt.xlim(1,T)
plt.legend()

plt.save("OscillatorError0,01.png")
