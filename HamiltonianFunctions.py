import numpy as np

def HaOscillator(tab_y):
    tab_p = tab_y[0]
    tab_q = tab_y[1]
    return 0.5*(np.power(tab_p,2) + np.power(tab_q,2))

def Oscillator(tn,yn):
    return np.array([-yn[1],yn[0]])

def Kepler(t,y):
    p = y[0:2]
    q = y[2:4]
    const = (q[0]**2 +q[1]**2)**(-3/2)
    return np.reshape(np.array([-q*const,p]),4)

def HaKepler(tab_y):
    tab_p1 = tab_y[:,0]
    tab_p2 = tab_y[:,1]
    tab_q1 = tab_y[:,2]
    tab_q2 = tab_y[:,3]
    return 0.5*(np.power(tab_p1,2) +np.power(tab_p2,2)) - np.power(np.power(tab_q1,2) +np.power(tab_q2,2),-1/2)

def Vlasov(t,y):
    x = y[:len(y)//2]
    v = y[len(y)//2:]
    E = np.array([-x[0],x[1]/2,x[2]/2])
    temp = np.array([0,v[2],-v[1]])
    return np.concatenate((np.copy(v),E+temp))   