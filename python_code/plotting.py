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

###############################
### plot the solutions
############################

def plot_Fokker_sphere(pi,level_set=None):
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

    #ax.plot_surface(x_sphere, y_sphere, z_sphere,  rstride=1, cstride=1, lw=0,cmap=plt.cm.Blues,alpha=0.5)

    #plot the solutions
    R = np.sqrt(pi[0][0,0]**2+pi[0][0,1]**2+pi[0][0,2]**2)
    for i in range(n_sol):
	#ax.plot(pi[i][:,0]/R,- pi[i][:,1]/R, pi[i][:,2]/R,c='k',lw=0.7)
	ax.scatter(pi[i][-1,0]/R,- pi[i][-1,1]/R, pi[i][-1,2]/R,c='k')
	ax.scatter(pi[i][0,0]/R,- pi[i][0,1]/R, pi[i][0,2]/R,c='r')

    plt.axis([-c,c,-c,c])
    ax.set_zlim3d(-c, c)
    ax.set_zlabel("$\Pi_3$")
    ax.set_ylabel("$\Pi_2$")
    ax.set_xlabel("$\Pi_1$")

def plot_distribution(distr,i):
    plt.figure()
    plt.imshow(np.transpose(distr),interpolation='none')
    plt.xlabel(r"$\phi$")
    plt.ylabel(r"$\theta$")
    plt.colorbar()
    plt.savefig("images/im_{0:s}".format(str(i).zfill(3)))
    plt.close()

def plot_momentum_sphere(pi,level_set=None):
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
    if level_set:
        piLS = level_sets()
        for i in range(len(piLS)):
	    R = np.sqrt(piLS[i][0,0]**2+piLS[i][0,1]**2+piLS[i][0,2]**2)
	    ax.plot(piLS[i][:,0]/R,- piLS[i][:,1]/R, piLS[i][:,2]/R,c='k',lw=0.5)

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


def plot_energy(K,C):
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

def plot_momentum(pi):
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

