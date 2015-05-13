import sys as sys
import numpy as np



###############################
### rigid body parameters  ####
###############################

I	= np.array([10.,20,30.]) #inertia tensor
dt      = 0.1    		#time resolution
t_max  = 10   #max time

T = np.arange(0,t_max,dt) #time vector, for very long run, we should avoid that, it uses to much memory


sig = 0.05 #noise width

n_sol = 200    #number of solutions

#fokker planck resolution

fokker=True
Nphi = 100
Ntheta = 50
