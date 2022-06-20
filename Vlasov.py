import numpy  as np


#Functions
def R(theta):
	mat = np.zeros((3,3))
	mat[0,0]=1
	mat[1,1] = np.cos(theta)
	mat[1,2] = np.sin(theta)
	mat[2,1] = -np.sin(theta)
	mat[2,2] = np.cos(theta)
	return mat

def Rcal(theta):
	mat = np.zeros((3,3))
	mat[1,1] = np.sin(theta)
	mat[1,2] = 1-np.cos(theta)
	mat[2,1] = np.cos(theta)-1
	mat[2,2] = np.sin(theta)
	return mat

def sol_ex(t,s,x,v,c,eps):
	gamma = eps/np.sqrt(1-2*c*eps**2)
	a = (1+np.sqrt(1+2*c*eps**2))/(2*eps)
	b = (1+np.sqrt(1-2*c*eps**2))/(2*eps)
	a1 = (-b*x[2]+v[1])*gamma
	a2 = (b*x[1]+v[2])*gamma
	b1 = (a*x[2]+v[1])*gamma
	b2 = (-a*x[1]-v[2])*gamma	
	c1 = x[0]
	c2 = v[0]/np.sqrt(c)
	ans = np.zeros(6)
	ans[0] = c1 * np.cos(np.sqrt(c) * (t-s)) + c2 * np.sin(np.sqrt(c)*(t-s))
	ans[1] = a1 * np.sin(a*(t-s)) - a2 * np.cos(a*(t-s)) + b1 * np.sin(b*(t-s)) - b2 * np.cos(b*(t-s))
	ans[2] = a1 * np.cos(a*(t-s)) + a2 * np.sin(a*(t-s)) + b1 * np.cos(b*(t-s)) + b2 * np.sin(b*(t-s))

	ans[3] = np.sqrt(c)*(c1 * -np.sin(np.sqrt(c) * (t-s)) + c2 * np.cos(np.sqrt(c)*(t-s)))
	ans[4] = a*(a1 * np.cos(a*(t-s)) + a2 * np.sin(a*(t-s))) + b*(b1 * np.cos(b*(t-s)) + b2 * np.sin(b*(t-s)))
	ans[5] = a*(-a1 * np.sin(a*(t-s)) + a2 * np.cos(a*(t-s))) + b*(-b1 * np.sin(b*(t-s)) + b2 * np.cos(b*(t-s)))

	return ans

def sol_ex_tab(tab_t,s,x,v,c,eps):
	tab_s = s*np.ones(len(tab_t))
	gamma = eps/np.sqrt(1-2*c*eps**2)
	a = (1+np.sqrt(1+2*c*eps**2))/(2*eps)
	b = (1+np.sqrt(1-2*c*eps**2))/(2*eps)
	a1 = (-b*x[2]+v[1])*gamma
	a2 = (b*x[1]+v[2])*gamma
	b1 = (a*x[2]+v[1])*gamma
	b2 = (-a*x[1]-v[2])*gamma	
	c1 = x[0]
	c2 = v[0]/np.sqrt(c)
	ans = np.zeros((len(tab_t),6))
	ans[:,0] = c1 * np.cos(np.sqrt(c) * (tab_t-tab_s)) + c2 * np.sin(np.sqrt(c)*(tab_t-tab_s))
	ans[:,1] = a1 * np.sin(a*(tab_t-tab_s)) - a2 * np.cos(a*(tab_t-tab_s)) + b1 * np.sin(b*(tab_t-tab_s)) - b2 * np.cos(b*(tab_t-tab_s))
	ans[:,2] = a1 * np.cos(a*(tab_t-tab_s)) + a2 * np.sin(a*(tab_t-tab_s)) + b1 * np.cos(b*(tab_t-tab_s)) + b2 * np.sin(b*(tab_t-tab_s))

	ans[:,3] = np.sqrt(c)*(c1 * -np.sin(np.sqrt(c) * (tab_t-tab_s)) + c2 * np.cos(np.sqrt(c)*(tab_t-tab_s)))
	ans[:,4] = a*(a1 * np.cos(a*(tab_t-tab_s)) + a2 * np.sin(a*(tab_t-tab_s))) + b*(b1 * np.cos(b*(tab_t-tab_s)) + b2 * np.sin(b*(tab_t-tab_s)))
	ans[:,5] = a*(-a1 * np.sin(a*(tab_t-tab_s)) + a2 * np.cos(a*(tab_t-tab_s))) + b*(-b1 * np.sin(b*(tab_t-tab_s)) + b2 * np.cos(b*(tab_t-tab_s)))
	return ans


def parareal_vlasov(F,g,x,v,s,tab_t,c,eps,delta_t,N,dtg):
    kmax = N
    y0 = np.concatenate((x,v))
    y_tab = np.zeros((kmax+2,N+1,len(y0)))
    lambda_f = np.zeros((len(y0),N+1))
    lambda_g = np.zeros((len(y0),N+1))
    y_tab[0] = sol_ex_tab(tab_t,s,x,v,c,eps) 
    for n in range(1,N+1):
        y_tab[1,n] = g(F, (n-1)*delta_t, n*delta_t, y_tab[1,n-1], dtg)
    k = 2
    while k<= kmax+1:
        y_tab[k,0] = y0
        for n in range(1,N+1):
            lambda_g[:,n]= g(F, (n-1)*delta_t, n*delta_t, y_tab[k-1,n-1], dtg)
        for n in range(1,N+1):
            y_tab[k,n] = lambda_f[:,n]-lambda_g[:,n] + g(F, (n-1)*delta_t, n*delta_t, y_tab[k,n-1], dtg)
        print("k:",k-1)    
        k+=1
    return y_tab    