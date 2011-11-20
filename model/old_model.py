#!/usr/bin/env pyth*Gon
#code design to calculate temperature and density profiles for PSR models of IP


import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from math import *
#from scipy.optimize import curve_fit
import stat

zzstep = 0.0112  # step je v Rwd !!!!!!

#const
while zzstep < 0.2 : 
	Rwd = 0.0112                            # polomer WD
	zzstep = zzstep + 0.001                  # tu sa nastavuje skript 
	z0 = 0.232				#toto je zle ked som tam dal 0 tak to zacalo hadzat velmi vysoke hodnoty 

	Mwd = 0.7                               # hmotnost WD 
	a = 1                                   # local mass acreation unit
	G = 6.6742e-11                          # gravitacna konstanta  
	mH = 1.66053886e-27                     # atomic mass const
	zstep = 0.00001                       # krok 
	u = 0.65                                # mean molecular weight  

	Msl =1.9891e30                          # hmotnost slnka  
	Rsl = 1.392e9                           # polomer slnka
	k = 1.3806505e-23		        # bolzmanova konstatna 

	K1 = (3*u*mH) / 1.3806505e-23           # konstanta v T0


	z = Rwd + zstep                         # pociatocne z
	count = 0
	dist = 0
	D = [] 
	T = []
	R = []

	v0 = 0.25 * sqrt((2*G*Mwd*Msl)/((z0*Rwd)*Rsl))                    # rychlost castice v zavislosti na z0 tu je chyba !!!
	T0 = K1*(v0**2)                                   # teplota shockovej vlny
	R0 = a / v0

	while dist < 0.015 : 
		
		T.append(  T0*(((z - Rwd)/((z0-Rwd)*Rwd))**0.4) )
		R.append(  R0*(((z - Rwd)/((z0-Rwd)*Rwd))**-0.4) )
		dist = (z - Rwd) / Rwd
		D.append(dist)
		count = count + 1
		z = z + zstep





#	print v0
	ttt = max(T) / 11604505
#	print T0
#	print K1
#	print max(T) 
	print z0


	#pocitanie spektra

#E = 5000*1.60217646e-19
	XX = []
	YY = []
	jz = []
	for ii in range (1000,500000):
		XX.append(ii)
		E = ii*1.60217646e-19
		jz = []
		for i in range (0,len(T)):
			j = 9.52e-38 * ((R[i] / (u*mH))**2) * (T[i]**-0.5) * ((E / (k*T[i]))**-0.4) * exp(-(E/(k*T[i])))	
			jz.append(j)	
		YY.append(sum(jz))



# normovanie kvoli plotu 
	Max = max(T)
	MaxR = max(R)
	for i in range(0,len(T)):
        	T[i] = T[i] / Max
        	R[i] = R[i] / MaxR

	fig = plt.figure(figsize=(8,8))
	ax = fig.add_subplot(311)
	fig.subplots_adjust(wspace=0.25, hspace=0.25, top=0.96, bottom=0.05, left=0.1, right=0.96)
	ax.plot(D,T,'-')
	ax.plot(D,R,'-')
	ax.grid(True)
	ax.set_yscale('log')
	ax.set_xlim((-0.002,0.017))
	ax.set_ylim((0,1.1))
	ax.set_xlabel("(z-Rwd)/Rwd")
	ax.set_ylabel("T/T_max , Ro/Ro_max")
	ax.legend(("T/T_max","R/R_max"),'upper right')


	bx = fig.add_subplot(312)
	bx.plot(XX,YY,"-")
	bx.grid(True)
	bx.set_yscale('log')
	bx.set_xscale('log')
	bx.set_xlabel("Energy (keV)")
	bx.set_ylabel("Flux")
	bx.set_xlim((10e2,10e5))	
	bx.set_ylim((1e-14,1e6))

	cx = fig.add_subplot(313)
	cx.text(0.1,0.1,'z0 = %4.4f Rwd' % (z0))
	cx.text(0.1,0.2,'v0 = %4.4f ms^-2' % (v0))
	cx.text(0.1,0.3,'T0 = %4.4f keV ' % (ttt))
	cx.set_xticks([])
        cx.set_yticks([])
	
	fpath = str(z0)

	plt.show()
#	plt.savefig(fpath,format='png', dpi=80,facecolor='w', edgecolor='w', orientation='portrait', papertype=None, transparent=False)
#	plt.close()
