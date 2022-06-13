# Import
import numpy as np
import matplotlib.pyplot as plt

# Parareal Algorithm


def parareal(F, f, g, y0, dtf, dtg, delta_t,T):
    N = int(T/delta_t)
    kmax = 10
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
        
        lambda_f[:,0] = np.copy(y0)
        lambda_tab[:, 0] = np.copy(y0)
        for n in range(1,N+1):
            lambda_f[:, n] = f(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtf) 
            lambda_g[:,n] = g(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtg)
        for n in range(1, N+1):
            lambda_tab[:, n] = np.copy(lambda_f[:,n])- np.copy(lambda_g[:,n])+  g(F, (n-1)*delta_t, n*delta_t, lambda_tab[:, n-1], dtg)
        lambda_tab_old = np.copy(lambda_tab)
        print("k:", k)
        np.linalg.norm(lambda_tab)
        k += 1
    return lambda_tab


def parareal_bis(F,f,g,y0,dtf,dtg,delta_t,T):
    kmax = 6
    N = int(T/delta_t)
    print("T:",T)
    print("N:",N)
    print("delta_t:",delta_t)
    print("dtg",dtg)
    print("dtf",dtf)
    y_tab = np.zeros((kmax+2,N+1,len(y0)))
    lambda_f = np.zeros((len(y0),N+1))
    
    y_tab[0,0] = y0 
    y_tab[1,0] = y0
    for n in range(1,N+1):
        y_tab[0,n] = f(F, (n-1)*delta_t, n*delta_t, y_tab[0,n-1], dtf)
        y_tab[1,n] = g(F, (n-1)*delta_t, n*delta_t, y_tab[1,n-1], dtg)
    k = 2
    while k<= kmax+1:
        y_tab[k,0] = y0
        for n in range(1,N+1):
            lambda_f[:,n] = f(F, (n-1)*delta_t, n*delta_t, y_tab[k-1,n-1], dtf)
            lambda_f[:,n]-= g(F, (n-1)*delta_t, n*delta_t, y_tab[k-1,n-1], dtg)
        for n in range(1,N+1):
            y_tab[k,n] = lambda_f[:,n] + g(F, (n-1)*delta_t, n*delta_t, y_tab[k,n-1], dtg)
        print("k:",k-1)    
        k+=1
    return y_tab    


def err(tab_y, sol_ex):
    K = len(tab_y)
    N = len(tab_y[0])
    tab_err = np.zeros((K,N))
    for k in range(K):
        for n in range(N):
            tab_err[k,n] = np.linalg.norm(tab_y[k,n]-sol_ex[:,n])
    return tab_err


def Ha_err(tab_y, tab_Ha0,Ha):
    K = len(tab_y)
    N = len(tab_y[0])
    tab_Ha = np.zeros((K,N))
    for k in range(K):
        tab_Ha[k] = tab_Ha0-Ha(tab_y[k])
    return np.abs(tab_Ha)


def parareal_Ha_Oscillator(F, f, g, T, y0, dtf, dtg, N, epsilon):
    delta_t = T/(N+1)
    kmax = 6
    lambda_tab_old = np.zeros((len(y0), N+1))
    lambda_f = np.zeros((len(y0),N+1))
    lambda_g = np.zeros((len(y0),N+1))
    lambda_tab = np.zeros((len(y0), N+1))
    tab_Ha = 0.5 * np.ones((kmax+2, N+1))
    tab_f = np.zeros((len(y0),N+1))

    tab_f[:,0] = np.copy(y0)

    for n in range(1,N+1):
        tab_f[:, n] = f(F, (n-1)*delta_t, n*delta_t, tab_f[:, n-1], dtf)
    tab_Ha[0] -= 0.5 *(np.power(tab_f[0],2) + np.power(tab_f[1],2))
    
    lambda_tab_old[:, 0] = np.copy(y0)
    for n in range(1, N+1):
        lambda_tab_old[:, n] = g(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtg)
    tab_Ha[1] -= 0.5 *(np.power(lambda_tab_old[0],2) + np.power(lambda_tab_old[1],2))
    
    k = 1
    lambda_tab[:, 0] = np.copy(y0)
    lambda_f[:,0] = np.copy(y0)
    lambda_g[:,0] = np.copy(y0)
    while k <= kmax:
        for n in range(1,N+1):
            lambda_f[:, n] = f(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtf) 
            lambda_g[:,n] = g(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtg)
        
        for n in range(1, N+1):
            lambda_tab[:,n] = lambda_f[:,n] - lambda_g[:,n] + g(F, (n-1)*delta_t, n*delta_t, lambda_tab[:, n-1], dtg)
        
        tab_Ha[k+1] -= 0.5 *(np.power(lambda_tab[0],2) + np.power(lambda_tab[1],2))

        lambda_tab_old = np.copy(lambda_tab)

        
        print("k:", k-1)
        k += 1
    return np.abs(tab_Ha)

def parareal_Ha_Kepler(F, f, g, T, y0, dtf, dtg, N, epsilon):
    delta_t = T/(N+1)
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