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
            lambda_f[:,n] = f(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtf) 
            lambda_g[:,n] = g(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtg)
        for n in range(1, N+1):
            lambda_tab[:, n] = np.copy(lambda_f[:,n])- np.copy(lambda_g[:,n])+  g(F, (n-1)*delta_t, n*delta_t, lambda_tab[:, n-1], dtg)
        lambda_tab_old = np.copy(lambda_tab)
        print("k:", k)
        np.linalg.norm(lambda_tab)
        k += 1
    return lambda_tab


def parareal_bis(F,f,g,y0,dtf,dtg,delta_t,T):
    kmax = 16
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


def parareal_bis_vlasov(F,f,g,y0,c,eps,delta_t,T,kmax):
    N = int(T/delta_t)
    print("T:",T)
    print("N:",N)
    print("delta_t:",delta_t)
    print("c:",c)
    print("eps:",eps)
    print("kmax:", kmax)
    y_tab = np.zeros((kmax+2,N+1,len(y0)))
    lambda_f = np.zeros((len(y0),N+1))
    y_tab[0,0] = y0 
    y_tab[1,0] = y0
    for n in range(1,N+1):
        y_tab[0,n] = f(F, (n-1)*delta_t, n*delta_t, y_tab[0,n-1], c,eps)
        y_tab[1,n] = g(F, (n-1)*delta_t, n*delta_t, y_tab[1,n-1], c,eps)
    k = 2
    while k<= kmax+1:
        y_tab[k,0] = y0
        for n in range(1,N+1):
            lambda_f[:,n] = f(F, (n-1)*delta_t, n*delta_t, y_tab[k-1,n-1], c,eps)
            lambda_f[:,n]-= g(F, (n-1)*delta_t, n*delta_t, y_tab[k-1,n-1], c,eps)
        for n in range(1,N+1):
            y_tab[k,n] = lambda_f[:,n] + g(F, (n-1)*delta_t, n*delta_t, y_tab[k,n-1], c, eps)
        print("k:",k-1)    
        k+=1
    return y_tab  


def parareal_bis_magnetic(F,f,g,y0,c,eps,delta_t,T,kmax):
    N = round(T/delta_t)
    print("T:",T)
    print("N:",N)
    print("delta_t:",delta_t)
    print("c:",c)
    print("eps:",eps)
    print("kmax:", kmax)
    y_tab = np.zeros((kmax+2,N+1,len(y0)))
    lambda_f = np.zeros((len(y0),N+1))
    y_tab[0,0] = y0 
    y_tab[1,0] = y0
    for n in range(1,N+1):
        y_tab[0,n] = f(F, (n-1)*delta_t, n*delta_t, y_tab[0,n-1], c,eps)
        y_tab[1,n] = g(F, (n-1)*delta_t, n*delta_t, y_tab[1,n-1], c,eps)
    k = 2
    while k<= kmax+1:
        y_tab[k,0] = y0
        for n in range(1,N+1):
            lambda_f[:,n] = f(F, (n-1)*delta_t, n*delta_t, y_tab[k-1,n-1], c,eps)
            lambda_f[:,n]-= g(F, (n-1)*delta_t, n*delta_t, y_tab[k-1,n-1], c,eps)
        for n in range(1,N+1):
            y_tab[k,n] = lambda_f[:,n] + g(F, (n-1)*delta_t, n*delta_t, y_tab[k,n-1], c, eps)
        print("k:",k-1)    
        k+=1
    return y_tab  

def parareal_Ha(F, f, g, y0, dtf, dtg, delta_t,T,Ha, tab_Ha0,tab_f):
    """Computes the hamiltonian of the system (works for Oscillator and Kepler)"""
    N = int(T/delta_t)
    kmax = 6
    lambda_tab_old = np.zeros((len(y0), N+1))
    lambda_tab = np.zeros((len(y0), N+1))
    lambda_f = np.zeros((len(y0),N+1))
    lambda_g = np.zeros((len(y0),N+1))
    tab_Ha = np.zeros((kmax+2,N+1))
    temp = tab_f(F,0,T,y0,delta_t)
    tab_Ha[0] = 0.5*(np.power(temp[0],2) + np.power(temp[1],2))-tab_Ha0
    lambda_tab_old[:, 0] = y0
    for n in range(1, N+1):
        lambda_tab_old[:, n] = g(F, (n-1)*delta_t, n*delta_t, lambda_tab_old[:, n-1], dtg)

    tab_Ha[1] = 0.5*(np.power(lambda_tab_old[0],2) + np.power(lambda_tab_old[1],2))-tab_Ha0
    
    lambda_f[:,0] = y0
    lambda_tab[:, 0] = y0
    lambda_g[:,0] = y0
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
        tab_Ha[k] = 0.5*(np.power(lambda_tab[0],2) + np.power(lambda_tab[1],2))-tab_Ha0
        print("k:", k)
        k += 1
    return np.abs(tab_Ha)

def err(tab_y, sol_ex):
    K = len(tab_y)
    N = len(tab_y[0])
    tab_err = np.zeros((K,N))
    for k in range(K):
        for n in range(N):
            tab_err[k,n] = np.linalg.norm(tab_y[k,n]-sol_ex[:,n])
    return tab_err


def Ha_err(tab_y, tab_Ha0,Ha):
    """Provided with the tab of Kth iterations of the parareal, computes the error
    in Hamiltonian for the system studied."""
    K = len(tab_y)
    N = len(tab_y[0])
    tab_Ha = np.zeros((K,N))
    for k in range(K):
        tab_Ha[k] = tab_Ha0-Ha(tab_y[k])
    return np.abs(tab_Ha)

