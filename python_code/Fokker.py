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
from utils  import *
from plotting  import *

###########################################
##compute solutions of the rigid body #####
###########################################

#initial conditions 

a_0 = 0.2*np.pi*abs(np.random.rand(n_sol,2))-[1,1]
#a_0 	= np.array([[2.,-1.],[2,-1],[2,-1]]) #initial momentum, in (\phi,\theta) coordinates

pi_0 = np.empty([n_sol,3]) 

for i in range(n_sol):
    pi_0[i]	= np.array([  np.sin(a_0[i,1]) * np.cos(a_0[i,0]),np.sin(a_0[i,1]) * np.sin(a_0[i,0]), np.cos(a_0[i,1]) ])

#integrate the solution
pi = []
for i in range(n_sol):
    pi.append(integration_heun_split(pi_0[i,:],T))

#compute the energy and casimirs

K = np.empty([n_sol,len(T)]) #energy
C = np.empty([n_sol,len(T)]) #casimir

for i in range(n_sol):
	#compute the kinetic energy
	K[i] = 0.5 * (pi[i][:,0]**2/I[0] + pi[i][:,1]**2/I[1] + pi[i][:,2]**2/I[2])
	#compute the casimir
	C[i] = 0.5 * (pi[i][:,0]**2 + pi[i][:,1]**2 + pi[i][:,2]**2)



plot_Fokker_sphere(pi)
for i in range(len(T)):
    plot_distribution(collect_fokker(pi,i),i)


#plot_momentum_sphere(pi)
#plot_energy(K,C)

plt.show()


