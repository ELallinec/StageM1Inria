##Imports
import numpy as np
import matplotlib.pyplot as plt


##Constants
t0 = 0
T= 10
y0 = np.array([1,0])

##Functions
def k1(f,tn,yn):
    return f(tn,yn)

def k2(f,tn,yn,h):
    return f(tn+h/2,yn+h*k1(f,tn,yn)/2)

def k3(f,tn,yn,h):
    return f(tn+h/2,yn+h*k2(f,tn,yn,h)/2)

def k4(f,tn,yn,h):
    return f(tn+h,yn+h*k3(f,tn,yn,h))

##f1 = J^-1 Grad(H(p,q)) avec H(p,q) = 1/2(p^2+q^2)
def f1(tn,yn):
    return np.array([-yn[1],yn[0]])

    

def RK4(f,t0,y0,h,m):
    y = y0
    tab_y = np.zeros((2,m))
    tab_Ha = np.zeros(m) 
    for k in range(m):
        tab_Ha[k] = (y[0]**2+y[1]**2)/2
        tab_y[:,k]=y
        y = y + (1/6) * h * ( k1(f,t0+h*k,y) + 2 * k2(f,t0+h*k,y,h) + 2 * k3(f,t0+h*k,y,h) + k4(f,t0+h*k,y,h) )   
    return tab_y, tab_Ha


#Error computation
tab_err = np.zeros((2,4))
tab_h = [0.1,0.01,0.001,0.0001]
k = 0
for h in tab_h:

    m = int(T/h)
    print(h,m,T)
    tab_t = np.linspace(t0,T,m)
    tab_y, Ha = RK4(f1,t0,y0,h,m)
    x_ex = np.cos(tab_t)
    v_ex = np.sin(tab_t)
    tab_err[:,k] = [np.max(np.abs(x_ex-tab_y[0])),np.max(np.abs(v_ex-tab_y[1]))]
    
    k+=1


print(tab_err)
plt.plot(np.log(tab_h),tab_err[0,:])
plt.plot(np.log(tab_h), tab_err[1,:])
plt.show()

