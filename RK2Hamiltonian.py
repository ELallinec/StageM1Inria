##Imports
import numpy as np
import matplotlib.pyplot as plt


##Constants
h = 0.01
m = 1000
y0 = np.array([1,0])

#Functions

def y_n(yn, tn,f,h):
    return yn +h*f(tn,yn)/2

def dy_n(yn,tn,f,h):
    return f(tn+h/2,y_n(yn,tn,f,h))

##f1 = J^-1 Grad(H(p,q)) avec H(p,q) = 1/2(p^2+q^2)
def f1(tn,yn):
    return np.array([-yn[1],yn[0]])

def RK2(f,t0,y0,h,m):
    y = y0
    tab_y = np.zeros((2,m))
    tab_Ha = np.zeros(m) 
    for k in range(m):
        tab_Ha[k] = (y[0]**2+y[1]**2)/2
        tab_y[:,k]=y
        y = y +h*dy_n(y,t0+k*h, f,h)
    return tab_y, tab_Ha


tab_y, tab_Ha = RK2(f1,0,y0,h,m)

#Plot approx vs exact solution
plt.plot(tab_y[0],tab_y[1])
plt.title("Phase portrait")
plt.xlabel("Position")
plt.ylabel("Velocity")
plt.show()

tab_t = np.linspace(0,m*h,m)
x_ex = np.cos(tab_t)
v_ex = np.sin(tab_t)

plt.plot(x_ex,v_ex)
plt.title("Exact solution")
plt.xlabel("Position") 
plt.ylabel("Velocity")
plt.show()
#Hamiltonian
plt.plot(tab_t,tab_Ha)
plt.title("Hamiltonian through time using RK2")
plt.show()


