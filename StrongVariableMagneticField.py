import numpy  as np
import RK4 as RK4

def SubSystem(t,y):
    y0 = y[:len(y)//2]
    u0 = y[len(y)//2:]
    A = np.zeros((3,3))
    A[0,0] = y0[1]**2
    A[0,1] = -y0[0]*y0[1]

    A[1,0] = -y0[0]*y0[1]
    A[1,1] = y0[0]**2

    A *= 1/(y0[0]**2+y0[1]**2)

    B = np.zeros(3)
    B[0] = u0[1]*(u0[0]*y0[1]-u0[1]*y0[0])/(y0[0]**2+y0[1]**2)
    B[1] = u0[0]*(u0[1]*y0[0]-u0[0]*y0[1])/(y0[0]**2+y0[1]**2)
    return np.concatenate((A.dot(u0),B))

def approx(s,t,y0,dt,eps):
    x = y0[:len(y0)//2]
    v = y0[len(y0)//2:]
    C = np.zeros((3,3))
    sol =  RK4.RK4(SubSystem,s,t,y0,dt)
    y = sol[:len(sol)//2]
    u = sol[len(sol)//2:]

    C[0, 0] = (y[0]**2*np.cos((t-s)/eps)+y[1]**2)/(y[0]**2+y[1]**2)
    C[0, 1] = y[0]*y[1]*(np.cos((t-s)/eps)-1)/(y[0]**2+y[1]**2)
    C[0, 2] = - y[0]*np.sin((t-s)/eps)/(np.sqrt(y[0]**2+y[1]**2))

    C[1, 0] = y[0]*y[1]*(np.cos((t-s)/eps)-1)/(y[0]**2+y[1]**2)
    C[1, 1] = (y[1]**2*np.cos((t-s)/eps)+y[0]**2)/(y[0]**2+y[1]**2)
    C[1, 2] = - y[1]*np.sin((t-s)/eps)/(np.sqrt(y[0]**2+y[1]**2))

    C[2, 0] = y[0]*np.sin((t-s)/eps)/(np.sqrt(y[0]**2+y[1]**2))
    C[2, 1] = y[1]*np.sin((t-s)/eps)/(np.sqrt(y[0]**2+y[1]**2))
    C[2, 2] = np.cos((t-s)/eps)

    G = np.concatenate((y,C.dot(u)))
    return G

def approx_tab(s,tab_t,y0,dt,eps):
    x = y0[:len(y0)//2]
    v = y0[len(y0)//2:]
    C = np.zeros((3,3))
    tab_y = np.zeros((len(tab_t),len(y0)))
    T = tab_t[-1]
    N = len(tab_t)-1
    n = 0   
    sol =  RK4.tab_RK4(SubSystem,s,T,y0,dt)
    for t in tab_t:
        y = sol[:len(sol)//2,n]
        u = sol[len(sol)//2:,n]
        C[0, 0] = (y[0]**2*np.cos((t-s)/eps)+y[1]**2)/(y[0]**2+y[1]**2)
        C[0, 1] = y[0]*y[1]*(np.cos((t-s)/eps)-1)/(y[0]**2+y[1]**2)
        C[0, 2] = - y[0]*np.sin((t-s)/eps)/(np.sqrt(y[0]**2+y[1]**2))

        C[1, 0] = y[0]*y[1]*(np.cos((t-s)/eps)-1)/(y[0]**2+y[1]**2)
        C[1, 1] = (y[1]**2*np.cos((t-s)/eps)+y[0]**2)/(y[0]**2+y[1]**2)
        C[1, 2] = - y[1]*np.sin((t-s)/eps)/(np.sqrt(y[0]**2+y[1]**2))

        C[2, 0] = y[0]*np.sin((t-s)/eps)/(np.sqrt(y[0]**2+y[1]**2))
        C[2, 1] = y[1]*np.sin((t-s)/eps)/(np.sqrt(y[0]**2+y[1]**2))
        C[2, 2] = np.cos((t-s)/eps)
        G = np.concatenate((y,C.dot(u)))
        tab_y[n] = G
        n +=1
    return tab_y


def parareal_bis_magnetic(F,f,g,y0,eps,dtf,delta_t,T,kmax):
    N = round(T/delta_t)
    print("T:",T)
    print("N:",N)
    print("delta_t:",delta_t)
    print("eps:",eps)
    print("kmax:", kmax)
    y_tab = np.zeros((kmax+2,N+1,len(y0)))
    lambda_f = np.zeros((len(y0),N+1))
    y_tab[0,0] = y0 
    y_tab[1,0] = y0
    for n in range(1,N+1):
        y_tab[0,n] = f(F, (n-1)*delta_t, n*delta_t, y_tab[0,n-1],dtf,eps)
        y_tab[1,n] = g((n-1)*delta_t, n*delta_t, y_tab[1,n-1], delta_t, eps)
    k = 2
    while k<= kmax+1:
        y_tab[k,0] = y0
        for n in range(1,N+1):
            lambda_f[:,n] = f(F, (n-1)*delta_t, n*delta_t, y_tab[k-1,n-1], dtf,eps)
            lambda_f[:,n]-= g((n-1)*delta_t, n*delta_t, y_tab[k-1,n-1], delta_t,eps)
        for n in range(1,N+1):
            y_tab[k,n] = lambda_f[:,n] + g( (n-1)*delta_t, n*delta_t, y_tab[k,n-1], delta_t, eps)
        print("k:",k-1)    
        k+=1
    return y_tab  

##RK4 using epsilon



def k1(f, tn, yn,eps):
    return f(tn, yn,eps)


def k2(f, tn, yn, dt,eps):
    return f(tn+dt/2, yn+dt*k1(f, tn, yn,eps)/2,eps)


def k3(f, tn, yn, dt,eps):
    return f(tn+dt/2, yn+dt*k2(f, tn, yn, dt,eps)/2,eps)


def k4(f, tn, yn, dt,eps):
    return f(tn+dt, yn+dt*k3(f, tn, yn, dt,eps),eps)

# f1 = J^-1 Grad(H(p,q)) avec H(p,q) = 1/2(p^2+q^2)


def small_RK4(f, t0, y0, dt,eps):
    return y0 + 1/6 * dt * (k1(f, t0+dt, y0,eps) + 2 * k2(f, t0+dt, y0, dt,eps) + 2 * k3(f, t0+dt, y0, dt,eps) + k4(f, t0+dt, y0, dt,eps))


def tab_RK4_magnetic(f, t0, T, y0, dt,eps):
    y = y0
    m = round((T-t0)/dt)
    tab_y = np.zeros((len(y0), m+1))
    tab_y[:, 0] = y0
    for k in range(m):
        y = small_RK4(f, t0+k*dt, y, dt,eps)
        tab_y[:, k+1] = y
    return tab_y

def tab_RK4_magnetic_2(f, t0, T, y0, dt,eps):
    y = y0
    m = round((T-t0)/dt)
    tab_y = np.zeros((len(y0), m+1))
    tab_y[:, 0] = y0
    for k in range(m):
        y = small_RK4(f, t0+k*dt-dt/2, y, dt,eps)
        tab_y[:, k+1] = y
    return tab_y

def RK4_magnetic(f, t0, T, y0, dt,eps):
    y = y0
    m = round((T-t0)/dt)
    for k in range(m):
        y = small_RK4(f, t0 + k*dt, y, dt,eps)
    return y



def G2(h,xk,vk,eps):
    bk = np.array([-xk[1]/np.sqrt(xk[0]**2+xk[1]**2),xk[0]/np.sqrt(xk[0]**2+xk[1]**2),eps])
    omega = -np.linalg.norm(bk,2)/eps 
    bk /= np.linalg.norm(bk,2)
    
    b = np.zeros((3,3))
    b[0,1] = -bk[2]
    b[0,2] = bk[1]
    
    b[1,0] = bk[2]
    b[1,2] = -bk[0]

    b[2,0] = -bk[1]
    b[2,1] = bk[0]
    E = np.eye(3) +np.sin(h*omega)*b + 0.5*(np.sin(h*omega/2)/0.5)**2 * (b@b)
    v =  E.dot(vk)
    x = xk + h*v
    return x,v


def G4(h,xk,vk,eps):
    g1 = 1/(2-2**(1/3))
    g0 = 1-2*g1
    x,v = xk,vk
    x,v = G2(g1*h,x,v,eps)
    x,v = G2(g0*h,x,v,eps)
    x,v = G2(g1*h,x,v,eps)
    return x,v

def G4_final(f,s,T,y0,dt,eps):
    x,v = y0[:len(y0)//2],y0[len(y0)//2:]
    m = round((T-s)/dt)
    for k in range(m):
        x,v = G4(dt,x,v,eps)
    return np.concatenate((x,v))

def G4_final_tab(f,s,T,y0,dt,eps):
    x,v = y0[:len(y0)//2],y0[len(y0)//2:]
    m = round((T-s)/dt)
    tab_y = np.zeros((6,m+1))
    tab_y[:,0] = np.concatenate((x,v))
    for k in range(m):
        tempv = v
        x,v = G4(dt,x,v,eps)
        tab_y[:,k+1] = np.concatenate((x,v))
    return tab_y