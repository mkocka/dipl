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
	Rwd = 7.8e8 * ((((1.44/Mwd)**2/3) - ((Mwd/1.44)**2/3))**1/2) 	
	Rwd = Rwd / 10   # returning in meters 
	return Rwd


def brem(Rwd,Mwd):
############################# constants ####################
        G = 6.67384e-11       # m^3 kg^-1 s^-2 
        u = 0.65              # mean molecular weight of fully ionized matter of solar abundance 
        mH = 1.66053886e-27   # atomic mass const
        k = 1.3806488e-23     # JK^-1 Boltzman const
	e = 1.60217646e-19
############################# magic constants #############
        a = 0.00001               # local mass acreation rate per area g cm^-2 s^-1 free param of model!!!
        z0 = 1.013          # Rwd units !!!! 
	zstep = 0.00001     # Rwd units !!!!
	Msun = 1.9891e30    # Sun mass kg 
###########################################################

	
	v0 = 0.25*sqrt((2*G*Mwd*Msun)/((z0)*Rwd))
	t0 = 3*((u*mH)/k)*(v0**2)
	rho0 = a/v0
	
	print v0,t0,rho0
	
	z = 1.00000001
	dist = []
	Ta = []
	Rhoa = []
	
		
 	while z < (z0) : 
		T = t0*(((z-1)/(z0-1))**0.4)  
		Rho = rho0*(((z-1)/(z0-1))**-0.4)
		dist.append(z-1)
		Ta.append(T)
		Rhoa.append(Rhoa)
		z = z + zstep	
		
	F = []	
	X= []
	Y = []
	for i in range(100,100000):
		E = i*e
		jz = []
		for ii in range(0,len(Ta)):
			j = 9.52e-38 * ((Rhoa[ii]/(u*mH))**2)*(Ta[ii]**-0.5)*((E/(k*Ta[ii]))**-0.4)*exp(-E/(k*Ta[ii]))
			jz.appned(j) 		
		X.append(E)
		Y.append(sum(jz))		

		
	return(Ta,Rhoa,dist,X,Y)

if __name__ == "__main__":
	if (len(sys.argv) == 2):
		Mwd = float(sys.argv[1])

		
	Rwd = Rwd_calc(Mwd)
	print Mwd,"Msun",Rwd,"m"

	ress = brem(Rwd,Mwd)
	#Ta,Rhoa,dist,X,Y
	






