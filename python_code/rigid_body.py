import sys as sys
import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator


from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif',size=15)

###############################
### rigid body parameters  ####
###############################

I	= np.array([10.,20,30.]) #inertia tensor
dt      = 0.1    		#time resolution
t_max  = 20000 #max time

sig = 0.01 #noise width

###############################
####numerical scheme##########3

def momentum(pi,t,i):
	pi_dot = np.empty(3)

	pi_dot[0] = pi[2]*pi[1]* (I[1]-I[2])/ (I[1]*I[2]) 
	pi_dot[1] = pi[0]*pi[2]* (I[2]-I[0])/ (I[2]*I[0])
	pi_dot[2] = pi[0]*pi[1]* (I[0]-I[1])/ (I[0]*I[1]) 
	return pi_dot 


def momentumS(pi,t,i):
	pi_dot = np.empty(3)

	pi_dot[0] = pi[2]*pi[1]* (I[1]-I[2])/ (I[1]*I[2]) 
	pi_dot[1] = pi[0]*pi[2]* (I[2]-I[0])/ (I[2]*I[0])
	pi_dot[2] = pi[0]*pi[1]* (I[0]-I[1])/ (I[0]*I[1]) 
	return pi_dot - sig**2*pi


#explicit Euler scheme
def integration_euler(pi_0, t):
    pi = np.empty([len(t),3])
    pi[0]=pi_0
    for i in range(len(t)-1):
        pi[i+1] = momentum(pi[i], t[i],i)*dt + pi[i]
    return pi

def integration_odeint(pi_0, t):
    pi = np.empty([len(t),3])
    pi[0]=pi_0


    for i in range(len(t)-1):
        pi[i+1] = integrate.odeint(momentum,pi[i,:],[0,dt],args=(i,))[1]
    return pi


# stochastic Euler 
def integration_eulerS(pi_0, t):
    pi = np.empty([len(t),3])
    pi[0]=pi_0


    np.random.seed(1)
    #dW=sig*np.sqrt(dt)*np.random.randn(3,len(T)) #random numbers
    dW=np.random.normal(0,sig*np.sqrt(dt),(3,len(T))) #random numbers

    for i in range(len(t)-1):
        pi[i+1] = momentumS(pi[i], t[i],i)*dt + pi[i] + np.cross(pi[i] ,dW[:,i]) -0.00*pi[i]*dt
    return pi


# stochastic split step
def integration_split(pi_0, t):
    pi = np.empty([len(t),3])
    pi[0]=pi_0


    np.random.seed(1)
    #dW=sig*np.sqrt(dt/1.)*np.random.randn(3,len(T)) #random numbers
    dW=np.random.normal(0,sig*np.sqrt(dt),(3,len(T))) #random numbers
    #dW[1,:]=0
    #dW[2,:]=0
    
    for i in range(len(t)-1):

        pi_tmp = integrate.odeint(momentumS,pi[i,:],[0,dt],args=(i,))[1]
        A = np.cross(pi_tmp ,dW[:,i])
        pi[i+1] =  pi_tmp +A# / np.linalg.norm(A)
        #pi[i+1] = pi[i+1]/np.linalg.norm(pi[i+1])
    return pi

# stochastic Heun for stratanovitch 
def integration_heun(pi_0, t):
    pi = np.empty([len(t),3])
    pi[0]=pi_0


    np.random.seed(1)
    #dW=sig*np.sqrt(dt/1.)*np.random.randn(3,len(T)) #random numbers
    dW=np.random.normal(0,sig*np.sqrt(dt),(3,len(T))) #random numbers
    for i in range(len(t)-1):
    	nu = pi[i] + dt*momentum(pi[i],t[i],i) + np.cross(pi[i] ,dW[:,i])
	pi[i+1] = pi[i] + dt*momentum(pi[i],t[i],i) + 0.5*( np.cross(pi[i] ,dW[:,i])+ np.cross(nu ,dW[:,i]) )
    return pi 

# stochastic Heun-split for stratanovitch 
def integration_heun_split(pi_0, t):
    pi = np.empty([len(t),3])
    pi[0]=pi_0


    np.random.seed(1)
    #dW=sig*np.sqrt(dt/1.)*np.random.randn(3,len(T)) #random numbers
    dW=np.random.normal(0,sig*np.sqrt(dt),(3,len(T))) #random numbers
    for i in range(len(t)-1):
        pi_tmp = integrate.odeint(momentum,pi[i,:],[0,dt],args=(i,))[1]
    	nu = pi_tmp +  np.cross(pi_tmp ,dW[:,i])
	pi[i+1] = pi_tmp + 0.5*( np.cross(pi_tmp ,dW[:,i])+ np.cross(nu ,dW[:,i]) )
    return pi 
##########################################
### compute the level sets  (DO NOT TOUCH!!)
#########################################

def level_sets():
    orbits =18 
    t_max   = [1000., 300.,200.,150.,120.,300.,200.,150.,120.,400.,300.,220.,200.,1000.,400.,300.,220.,200.] #max time

    a_0 	= np.array([[0.001,np.pi/2.],[0.0,np.pi/2.-0.3],[0.0,np.pi/2.-0.6],[0.0,np.pi/2.-0.9],[0.0,np.pi/2.-1.2],[0.0,np.pi/2.+0.3],[0.0,np.pi/2.+0.6],[0.0,np.pi/2.+0.9],[0.0,np.pi/2.+1.2],[0.3,np.pi/2.],[0.6,np.pi/2.],[0.9,np.pi/2.],[1.2,np.pi/2.],[-0.001,np.pi/2.],[-0.3,np.pi/2.],[-0.6,np.pi/2.],[-0.9,np.pi/2.],[-1.2,np.pi/2.]]) #initial positions


    pi_0 = np.empty([orbits,3])
    t=[]
    for i in range(orbits):
	pi_0[i]	= np.array([  np.cos(a_0[i,0]) * np.cos(a_0[i,1]),np.cos(a_0[i,0]) * np.sin(a_0[i,1]), np.sin(a_0[i,0]) ])
	t.append(np.arange(0,t_max[i],dt))

    pi = []
    for i in range(orbits):
	#pi.append(integrate.odeint(momentum,pi_0[i,:],t[i],args=(i,)))
	pi.append(integration_euler(pi_0[i,:],t[i]))

    return pi



###########################################
##compute solutions of the rigid body #####
###########################################

n_sol = 3    #number of solutions

T = np.arange(0,t_max,dt) #time vector, for very long run, we should avoid that, it uses to much memory


#initial conditions 
a_0 	= np.array([[2.,-1.],[2,-1],[2,-1]]) #initial momentum, in (\phi,\theta) coordinates

pi_0 = np.empty([n_sol,3]) 

for i in range(n_sol):
    pi_0[i]	= np.array([  np.cos(a_0[i,0]) * np.cos(a_0[i,1]),np.cos(a_0[i,0]) * np.sin(a_0[i,1]), np.sin(a_0[i,0]) ])

#integrate the solution
pi = []
for i in range(n_sol-2):
    #pi.append(integrate.odeint(momentum,pi_0[i,:],T,args=(i,)))
    #pi.append(integration_eulerS(pi_0[i,:],T))
    pi.append(integration_split(pi_0[i,:],T))
    pi.append(integration_heun_split(pi_0[i+1,:],T))
    pi.append(integration_odeint(pi_0[i+2,:],T))

#compute the energy and casimirs

K = np.empty([n_sol,len(T)]) #energy
C = np.empty([n_sol,len(T)]) #casimir

for i in range(n_sol):
	#compute the kinetic energy
	K[i] = 0.5 * (pi[i][:,0]**2/I[0] + pi[i][:,1]**2/I[1] + pi[i][:,2]**2/I[2])
	#compute the casimir
	C[i] = 0.5 * (pi[i][:,0]**2 + pi[i][:,1]**2 + pi[i][:,2]**2)



###############################
### plot the solutions
############################

#create the momentum sphere
c = 1. #radius
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

x_sphere = c * np.outer(np.cos(u), np.sin(v))
y_sphere = c * np.outer(np.sin(u), np.sin(v))
z_sphere = c * np.outer(np.ones(np.size(u)), np.cos(v))


#plot the sphere
fig = plt.figure(figsize=(8,7))
ax = fig.gca(projection='3d')

ax.plot_surface(x_sphere, y_sphere, z_sphere,  rstride=1, cstride=1, lw=0,cmap=plt.cm.Blues,alpha=0.5)


#plot the level sets
"""
piLS = level_sets()
for i in range(len(piLS)):
	R = np.sqrt(piLS[i][0,0]**2+piLS[i][0,1]**2+piLS[i][0,2]**2)
	ax.plot(piLS[i][:,0]/R,- piLS[i][:,1]/R, piLS[i][:,2]/R,c='k',lw=0.5)
"""
#plot the solutions
for i in range(n_sol-1):
	R = np.sqrt(pi[i][0,0]**2+pi[i][0,1]**2+pi[i][0,2]**2)
	ax.plot(pi[i][:,0]/R,- pi[i][:,1]/R, pi[i][:,2]/R,c='r',lw=1)
	ax.plot(pi[i+1][:,0]/R,- pi[i+1][:,1]/R, pi[i+1][:,2]/R,c='b',lw=1)

plt.axis([-c,c,-c,c])
ax.set_zlim3d(-c, c)
ax.set_zlabel("$\Pi_3$")
ax.set_ylabel("$\Pi_2$")
ax.set_xlabel("$\Pi_1$")


#plot kinetic and casimir energy
plt.figure()
for i in range(n_sol-2):
	plt.semilogy(T,abs((K[i,:]-K[i,0])/K[i,0]),label="$\Delta E,\ \mathrm{Ito}$")
	plt.semilogy(T,abs((C[i,:]-C[i,0])/C[i,0]),label="$\Delta C,\ \mathrm{Ito}$")
	#plt.plot(T,(C[i,:]-C[i,0])/C[i,0],label="$C(t)-C(0)\, Ito$")
	plt.semilogy(T,abs((K[i+1,:]-K[i+1,0])/K[i+1,0]),label="$\Delta E,\ \mathrm{Str}$")
	plt.semilogy(T,abs((C[i+1,:]-C[i+1,0])/C[i+1,0]),label="$\Delta C,\ \mathrm{Str}$")
	plt.semilogy(T,abs((K[i+2,:]-K[i+1,0])/K[i+1,0]),label="$\Delta E,\ \mathrm{Det}$")
	plt.semilogy(T,abs((C[i+2,:]-C[i+1,0])/C[i+1,0]),label="$\Delta C,\ \mathrm{Det}$")

plt.xlabel(r"$\mathrm{time}$")
plt.ylabel(r"$\mathrm{error}$")
plt.legend(loc='lower right')
plt.savefig("energy.png")

plt.figure()
plt.plot(T,pi[0][:,0],'r',label=r'$\Pi_1,\ \mathrm{Ito}$')
plt.plot(T,pi[1][:,0],'r--',label=r'$\Pi_1,\ \mathrm{Str}$')
plt.plot(T,pi[2][:,0],'r-.',label=r'$\Pi_1,\ \mathrm{Det}$')
plt.plot(T,pi[0][:,1],'b',label=r'$\Pi_2,\ \mathrm{Ito}$')
plt.plot(T,pi[1][:,1],'b--',label=r'$\Pi_2,\ \mathrm{Str}$')
plt.plot(T,pi[2][:,1],'b-.',label=r'$\Pi_2,\ \mathrm{Det}$')
plt.plot(T,pi[0][:,2],'g',label=r'$\Pi_3, \mathrm{Ito}$')
plt.plot(T,pi[1][:,2],'g--',label=r'$\Pi_3,\ \mathrm{Str}$')
plt.plot(T,pi[2][:,2],'g-.',label=r'$\Pi_3,\ \mathrm{Det}$')

plt.xlabel(r"$\mathrm{time}$")
plt.ylabel(r"$\Pi$")
plt.legend(loc='best')
plt.savefig("momentum.png")
#plt.close()

plt.show()
