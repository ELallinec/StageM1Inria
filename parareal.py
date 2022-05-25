# Import
import numpy as np
import matplotlib.pyplot as plt

# Parareal Algorithm


def parareal(F, f, g, T, y0, dtf, dtg, N, epsilon):
    delta_t = T/N
    kmax = N
    lambda_tab_old = np.zeros((len(y0), N+1))
    lambda_tab = np.zeros((len(y0), N+1))
    lambda_f = np.zeros((len(y0),N+1))
    lambda_g = np.zeros((len(y0),N+1))
    
    lambda_tab_old[:, 0] = np.copy(y0)
    for n in range(1, N+1):
        lambda_tab_old[:, n] = g(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtg)
    
    

    lambda_f[:,0] = np.copy(y0)
    lambda_tab[:, 0] = np.copy(y0)
    lambda_g[:,0] = np.copy(y0)
    k = 1
    while k <= kmax:
        print("k:", k)
        lambda_f[:,0] = np.copy(y0)
        lambda_tab[:, 0] = np.copy(y0)
        for n in range(1,N+1):
            lambda_f[:, n] = f(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtf) 
            lambda_g[:,n] = g(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtg)
        print("f",lambda_f)
        print("g",lambda_g)
        for n in range(1, N+1):
            print(f"lambda_tab {n} début\n", lambda_tab)
            
            lambda_tab[:, n] = np.copy(lambda_f[:,n])- np.copy(lambda_g[:,n])+  g(F, (n-1)*delta_t, n*delta_t, lambda_tab[:, n-1], dtg)
            print(f"lambda_tab {n} fin\n", lambda_tab)
        lambda_tab_old = np.copy(lambda_tab)
        
        
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
    lambda_tab_old = np.zeros((len(y0), N+1),dtype=np.double)
    lambda_f = np.zeros((len(y0),N+1),dtype=np.double)
    lambda_g = np.zeros((len(y0),N+1),dtype=np.double)
    lambda_tab = np.zeros((len(y0), N+1),dtype=np.double)
    tab_Ha = 0.5 * np.ones((kmax+2, N+1),dtype=np.double)
    tab_f = np.zeros((len(y0),N+1),dtype=np.double)

    tab_f[:,0] = np.copy(y0)

    for n in range(1,N+1):
        tab_f[:, n] = f(F, (n-1)*delta_t, n*delta_t, tab_f[:, n-1], dtf)
    tab_Ha[0] = 0.5 *(np.power(tab_f[0],2) + np.power(tab_f[1],2))
    
    lambda_tab_old[:, 0] = np.copy(y0)
    for n in range(1, N+1):
        lambda_tab_old[:, n] = g(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtg)
    tab_Ha[1] = 0.5 *(np.power(lambda_tab_old[0],2) + np.power(lambda_tab_old[1],2))
    
    k = 1
    lambda_tab[:, 0] = np.copy(y0)
    lambda_f[:,0] = np.copy(y0)
    lambda_g[:,0] = np.copy(y0)
    while k <= kmax:
        for n in range(1,N+1):
            lambda_f[:, n] = f(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtf) 
            lambda_g[:,n] = g(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtg)
        
        for n in range(1, N+1):
            lambda_tab[:,n] = np.copy(lambda_f[:,n]) - np.copy(lambda_g[:,n]) + g(F, (n-1)*delta_t, n*delta_t, lambda_tab[:, n-1], dtg)
        
        tab_Ha[k+1] = 0.5 *(np.power(lambda_tab[0],2) + np.power(lambda_tab[1],2))

        lambda_tab_old = np.copy(lambda_tab)

        
        print("k:", k)
        k += 1
    return lambda_tab

def parareal_Ha_Kepler(F, f, g, T, y0, dtf, dtg, N, epsilon):
    delta_t = T/N
    kmax = 6
    lambda_tab_old = np.zeros((len(y0), N+1),dtype=np.double)
    lambda_f = np.zeros((len(y0),N+1),dtype=np.double)
    lambda_tab = np.zeros((len(y0), N+1),dtype=np.double)
    tab_Ha = 0.5 * np.ones((kmax+2, N+1),dtype=np.double)
    tab_f = np.zeros((len(y0),N+1),dtype=np.double)

    tab_f[:,0] = y0
    tab_Ha[0,0]= 0
    for n in range(1,N+1):
        tab_f[:, n] = f(F, (n-1)*delta_t, n*delta_t, tab_f[:, n-1], dtf)
    tab_Ha[0] = 0.5 *(np.power(tab_f[0],2) + np.power(tab_f[1],2))- 1/np.sqrt(np.power(tab_f[2],2)+np.power(tab_f[3],2))
    
    lambda_tab_old[:, 0] = y0
    tab_Ha[1,0] = 0
    for n in range(1, N+1):
        lambda_tab_old[:, n] = g(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtg)
    tab_Ha[1] = 0.5 *(np.power(lambda_tab_old[0],2) + np.power(lambda_tab_old[1],2))- 1/np.sqrt(np.power(lambda_tab_old[2],2)+np.power(lambda_tab_old[3],2))
    
    k = 1
    lambda_tab[:, 0] = y0
    lambda_f[:,0] = y0
    while k <= kmax:
        for n in range(1,N+1):
            lambda_f[:, n] = f(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtf) 
        for n in range(1, N+1):
            lambda_tab[:,n] = lambda_f[:,n] - g(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtg)
            lambda_tab[:, n] +=  g(F, (n-1)*delta_t, n*delta_t, lambda_tab[:, n-1], dtg)
        tab_Ha[k+1] = 0.5 *(np.power(lambda_tab[0],2) + np.power(lambda_tab[1],2))- 1/np.sqrt(np.power(lambda_tab[2],2)+np.power(lambda_tab[3],2))
        lambda_tab_old = np.copy(lambda_tab)
        
        print("k:", k)
        k += 1
    return tab_Ha