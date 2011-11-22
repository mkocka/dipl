#!/usr/bin/env python
#code design to calculate temperature and density profiles for PSR models of IP

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from math import *
import stat

def Rwd_calc(Mwd):
	#Msun = 1.9891e30
	Rwd = 7.8e8 * ((((1.44/Mwd)**(2.0/3)) - ((Mwd/1.44)**(2.0/3)))**(1.0/2)) 	
	Rwd = Rwd / 100   # returning in meters 
	return Rwd


def brem(Rwd,Mwd):
############################# constants ####################
        G = 6.67384e-11       # m^3 kg^-1 s^-2 
        u = 0.65              # mean molecular weight of fully ionized matter of solar abundance 
        mH = 1.66053886e-27   # atomic mass const
        k = 1.3806488e-23     # JK^-1 Boltzman const
	e = 1.60217646e-19
############################# magic constants #############
        a = 1               # local mass acreation rate per area g cm^-2 s^-1 free param of model!!!
        z0 = 1.013          # Rwd units !!!! 
	zstep = 0.00001     # Rwd units !!!!
	Msun = 1.9891e30    # Sun mass kg 
###########################################################

	
	v0 = 0.25*sqrt((2*G*Mwd*Msun)/((z0)*Rwd))
	t0 = 3*((u*mH)/k)*(v0**2)
	rho0 = a/v0
	
	z = 1.00000001
	dist = []
	Ta = []
	Rhoa = []
	
		
 	while z < (z0) : 
		T = t0*(((z-1)/(z0-1))**0.4)  
		Rho = rho0*(((z-1)/(z0-1))**-0.4)
		dist.append(z-1)
		Ta.append(T)
		Rhoa.append(Rho)
		z = z + zstep	

	print (max(Ta)/11604505),"keV", t0
#	plt.plot(dist,Ta)
#	plt.plot(dist,Rhoa)
#	plt.show()	

#	print (max(Ta)/11604505),"keV"
	F = []	
	X= []
	Y = []
	MM = 11604505 # keV magic const
	i = 100
	while i < 1000000:
		E = i/1000.0
	#	print i
		jz = []
		for ii in range(0,len(Ta)):
			j = 9.52e-38 * ( (Rhoa[ii]/(u*mH))**(2) ) * (Ta[ii]**(-0.5)) * ( (E/(k*Ta[ii]/e))**(-0.4)) * exp(-E/(k*Ta[ii]/e))

	#		c1 =  (Rhoa[ii]/(u*mH))**(2)  
	#		c2 = Ta[ii]**(-0.5)  
	#		c3 =  (E/(k*Ta[ii]/e))**(-0.4) 
	#		c4 =  exp(-E/(k*Ta[ii]/e))

	#		print c1,c2,c3,c4
			jz.append(j) 		
		X.append(E)
		Y.append(sum(jz))		
		i = i + 500
		
	return(Ta,Rhoa,dist,X,Y)

def plot_temp_dens(T,Rho,D,X,Y):
	T_max = max(T)
	Rho_max = max(Rho)
	for i in range(0,len(T)):
		T[i] = T[i] / T_max
		Rho[i]= Rho[i] / Rho_max
	
	fig = plt.figure()
	ax=fig.add_subplot(211)
	ax.plot(D,T,'-')
	ax.plot(D,Rho,'-')
	ax.grid(True)
	ax.set_yscale('log')
        ax.set_xlabel("(z-Rwd)/Rwd")
        ax.set_ylabel("T/T_max , Ro/Ro_max")
        ax.legend(("T/T_max","R/R_max"),'upper right')
        ax.set_xlim((-0.002,0.017))
        ax.set_ylim((0,1.1))
	
	bx = fig.add_subplot(212)
	bx.plot(X,Y,'-')
        bx.grid(True)
        bx.set_yscale('log')
        bx.set_xscale('log')
        bx.set_xlabel("Energy (keV)")
        bx.set_ylabel("Flux")


	plt.show()
	plt.close()

	return()
	
	


if __name__ == "__main__":
	if (len(sys.argv) == 2):
		Mwd = float(sys.argv[1])

		
	Rwd = Rwd_calc(Mwd)
	print Mwd,"Msun",Rwd,"m"
	ress = brem(Rwd,Mwd)
	#Ta,Rhoa,dist,X,Y
	plot_temp_dens(ress[0],ress[1],ress[2],ress[3],ress[4])	
	





