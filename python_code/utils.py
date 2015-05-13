import sys as sys
import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator


from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif',size=15)

from param import *
from plotting import *
###############################
####numerical schemes##########

#without dissipation 
def momentum(pi,t,i):
	pi_dot = np.empty(3)

	pi_dot[0] = pi[2]*pi[1]* (I[1]-I[2])/ (I[1]*I[2]) 
	pi_dot[1] = pi[0]*pi[2]* (I[2]-I[0])/ (I[2]*I[0])
	pi_dot[2] = pi[0]*pi[1]* (I[0]-I[1])/ (I[0]*I[1]) 
	return pi_dot 


#with dissipation
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
    dW=np.random.normal(0,sig*np.sqrt(dt),(3,len(T))) #random numbers

    for i in range(len(t)-1):
        pi[i+1] = momentumS(pi[i], t[i],i)*dt + pi[i] + np.cross(pi[i] ,dW[:,i]) -0.00*pi[i]*dt
    return pi


# stochastic split step
def integration_split(pi_0, t):
    pi = np.empty([len(t),3])
    pi[0]=pi_0


    np.random.seed(1)
    dW=np.random.normal(0,sig*np.sqrt(dt),(3,len(T))) #random numbers
    
    for i in range(len(t)-1):

        pi_tmp = integrate.odeint(momentumS,pi[i,:],[0,dt],args=(i,))[1]
        pi[i+1] =  pi_tmp+ np.cross(pi_tmp ,dW[:,i])
    return pi

# stochastic Heun for stratanovitch 
def integration_heun(pi_0, t):
    pi = np.empty([len(t),3])
    pi[0]=pi_0

    np.random.seed(1)
    dW=np.random.normal(0,sig*np.sqrt(dt),(3,len(T))) #random numbers

    for i in range(len(t)-1):
    	nu = pi[i] + dt*momentum(pi[i],t[i],i) + np.cross(pi[i] ,dW[:,i])
	pi[i+1] = pi[i] + dt*momentum(pi[i],t[i],i) + 0.5*( np.cross(pi[i] ,dW[:,i])+ np.cross(nu ,dW[:,i]) )
    return pi 

# stochastic Heun-split for stratanovitch 
def integration_heun_split(pi_0, t):
    pi = np.empty([len(t),3])
    pi[0]=pi_0


    #np.random.seed(1)
    dW=np.random.normal(0,sig*np.sqrt(dt),(3,len(T))) #random numbers

    for i in range(len(t)-1):
        pi_tmp = integrate.odeint(momentum,pi[i,:],[0,dt],args=(i,))[1]
    	nu = pi_tmp +  np.cross(pi_tmp ,dW[:,i])
	pi[i+1] = pi_tmp + 0.5*( np.cross(pi_tmp ,dW[:,i])+ np.cross(nu ,dW[:,i]) )

    return pi 


#################################
###collect data for fokker-planck
############################

def collect_fokker(pi,tt):

    dphi = 2*np.pi/Nphi
    dtheta = np.pi/Ntheta
    Phi = np.arange(0,2*np.pi,dphi)
    Theta = np.arange(0,np.pi,dtheta)

    distr =np.zeros([Nphi,Ntheta])

    for i in range(n_sol):
        r = np.sqrt(pi[i][tt,0]**2+pi[i][tt,1]**2+pi[i][tt,2]**2)
        phi = np.arctan2(pi[i][tt,1],pi[i][tt,0])+np.pi
        theta = np.arccos(pi[i][tt,2]/r)
        for k in range(Nphi-1):
            if phi< Phi[k+1] and phi>= Phi[k]:
                for l in range(Ntheta-1):
                    if theta< Theta[l+1] and theta>= Theta[l]:
                        distr[k,l] +=1

    return distr

###################################
### compute the level sets ######## 
###################################

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



