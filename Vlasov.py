import numpy  as np


#Functions
#functions
def R(theta):
	mat = np.zeros((3,3))
	mat[0,0] = 1
	mat[1,1] = np.cos(theta)
	mat[1,2] = np.sin(theta)
	mat[2,2] = np.cos(theta)
	mat[2,1] = -np.sin(theta)
	return mat

def Rcal(theta):
	mat = np.zeros((3,3))
	mat[1,1] = np.sin(theta)
	mat[1,2] = 1-np.cos(theta)
	mat[2,2] = np.sin(theta)
	mat[2,1] = np.cos(theta)-1
	return mat

def sol_ex(F,s,t,y0,c,eps):

	x = y0[:3]
	v = y0[3:]
	gamma = eps/np.sqrt(1-2*c*eps**2)
	a = (1+np.sqrt(1-2*c*eps**2))/(2*eps)
	b = (1-np.sqrt(1-2*c*eps**2))/(2*eps)
	a1 = (-b*x[2]+v[1])*gamma
	a2 = (b*x[1]+v[2])*gamma
	b1 = (a*x[2]-v[1])*gamma
	b2 = (-a*x[1]-v[2])*gamma	
	c1 = x[0]
	c2 = v[0]/np.sqrt(c)
	ans = np.zeros(6)
	ans[0] = c1 * np.cos(np.sqrt(c) * (t-s)) + c2 * np.sin(np.sqrt(c)*(t-s))
	ans[1] = a1 * np.sin(a*(t-s)) - a2 * np.cos(a*(t-s)) + b1 * np.sin(b*(t-s)) - b2 * np.cos(b*(t-s))
	ans[2] = a1 * np.cos(a*(t-s)) + a2 * np.sin(a*(t-s)) + b1 * np.cos(b*(t-s)) + b2 * np.sin(b*(t-s))

	ans[3] = np.sqrt(c)*(-c1 * np.sin(np.sqrt(c) * (t-s)) + c2 * np.cos(np.sqrt(c)*(t-s)))
	ans[4] = a*(a1 * np.cos(a*(t-s)) + a2 * np.sin(a*(t-s))) + b*(b1 * np.cos(b*(t-s)) + b2 * np.sin(b*(t-s)))
	ans[5] = a*(-a1 * np.sin(a*(t-s)) + a2 * np.cos(a*(t-s))) + b*(-b1 * np.sin(b*(t-s)) + b2 * np.cos(b*(t-s)))
	return ans


def sol_ex_mat(F,s,t,y0,c,eps):
	gamma = eps/np.sqrt(1-2*c*eps**2)
	a = (1+np.sqrt(1-2*c*eps**2))/(2*eps)
	b = (1-np.sqrt(1-2*c*eps**2))/(2*eps)
	matA = np.longdouble(np.zeros((6,6)))
	matA[0,0] = np.cos(np.sqrt(c)*(t-s))
	matA[1,1] = gamma*(a*np.cos(b*(t-s)) -b*np.cos(a*(t-s)))
	matA[1,2] =  gamma*(a*np.sin(b*(t-s)) -b*np.sin(a*(t-s)))
	matA[2,1] =  gamma*(b*np.sin(a*(t-s)) -a*np.sin(b*(t-s)))
	matA[2,2] =  gamma*(a*np.cos(b*(t-s)) -b*np.cos(a*(t-s)))
	
	matA[0,3] = np.sin(np.sqrt(c)*(t-s))/np.sqrt(c)
	matA[1,4] = gamma*(np.sin(a*(t-s)) -np.sin(b*(t-s)))
	matA[1,5] =  gamma*(np.cos(b*(t-s)) - np.cos(a*(t-s)))
	matA[2,4] =  gamma*(np.cos(a*(t-s)) - np.cos(b*(t-s)))
	matA[2,5] =  gamma*(np.sin(a*(t-s)) -np.sin(b*(t-s)))
	
	matA[3,0] = -np.sqrt(c)*np.sin(np.sqrt(c)*(t-s))
	matA[4,1] = gamma*(-a*b*np.sin(b*(t-s)) +a*b*np.sin(a*(t-s)))
	matA[4,2] =  gamma*(a*b*np.cos(b*(t-s)) -a*b*np.cos(a*(t-s)))
	matA[5,1] =  gamma*(a*b*np.cos(a*(t-s)) -a*b*np.cos(b*(t-s)))
	matA[5,2] =  gamma*(-a*b*np.sin(b*(t-s)) + a*b*np.sin(a*(t-s)))
	
	matA[3,3] = np.cos(np.sqrt(c)*(t-s))
	matA[4,4] = gamma*(a*np.cos(a*(t-s)) -b*np.cos(b*(t-s)))
	matA[4,5] =  gamma*(-b*np.sin(b*(t-s)) + a*np.sin(a*(t-s)))
	matA[5,4] =  gamma*(-a*np.sin(a*(t-s)) +b*np.sin(b*(t-s)))
	matA[5,5] =  gamma*(a*np.cos(a*(t-s)) -b*np.cos(b*(t-s)))
	
	y = matA.dot(y0)
	return y

def sol_approx(F,s,t,y0,c,eps):
	P = np.longdouble(np.zeros((3,3)))
	P[0,0] = 1
	matA = np.longdouble(np.zeros((6,6)))
	matA[:3,:3] = np.eye(3)+c*eps*(t-s)/2 * R(np.pi/2)+(np.cos(np.sqrt(c)*(t-s))-1-c*eps*(t-s)/2)*P
	matA[:3,3:] = eps * Rcal((t-s)/eps) + np.sin(np.sqrt(c)*(t-s))/np.sqrt(c) * P
	matA[3:,:3] = eps*c/2 * Rcal((t-s)/eps) - np.sin(np.sqrt(c)*(t-s))*np.sqrt(c)*P
	matA[3:,3:] = R((t-s)/eps)+c*eps*(t-s)/2 *R((t-s)/eps-np.pi/2)+(np.cos(np.sqrt(c)*(t-s))-1-c*eps*(t-s)/2)*P
	y = matA.dot(y0)
	return ybm



def sol_ex_tab(tab_t,s,x,v,c,eps):
	tab_s = s*np.ones(len(tab_t))
	gamma = eps/np.sqrt(1-2*c*eps**2)
	a = (1+np.sqrt(1-2*c*eps**2))/(2*eps)
	b = (1-np.sqrt(1-2*c*eps**2))/(2*eps)
	a1 = (-b*x[2]+v[1])*gamma
	a2 = (b*x[1]+v[2])*gamma
	b1 = (a*x[2]-v[1])*gamma
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

def sol_approx_tab(tab_t,s,x,v,c,eps):
	P = np.zeros((3,3))
	P[0,0] = 1
	matA = np.zeros((6,6))
	y0 =np.concatenate((x,v))
	y_tab = np.zeros((len(tab_t),len(y0)))
	n =0
	for t in tab_t:
		matA = np.zeros((6,6))	
		matA[:3,:3] = np.eye(3)+c*eps*(t-s)/2 * R(np.pi/2)+(np.cos(np.sqrt(c)*(t-s))-1-c*eps*(t-s)/2)*P
		matA[:3,3:] = eps * Rcal((t-s)/eps) + np.sin(np.sqrt(c)*(t-s))/np.sqrt(c) * P
		matA[3:,:3] = eps/2 * Rcal((t-s)/eps) - np.sin(np.sqrt(c)*(t-s))*np.sqrt(c)*P
		matA[3:,3:] = R((t-s)/eps)+c*eps*(t-s)/2 *R((t-s)/eps-np.pi/2)+(np.cos(np.sqrt(c)*(t-s))-1-c*eps*(t-s)/2)*P
		y_tab [n] = matA.dot(y0)
		n += 1
	return y_tab


