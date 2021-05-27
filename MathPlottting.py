#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 26 15:14:38 2021

@author: 
   ___  _ ____             __          
  / _ \(_) / /__ ___      / /  ___ ___ 
 / // / / / / -_) _ \    / /__/ -_) -_)
/____/_/_/_/\__/_//_/   /____/\__/\__/ 

"""
#import the packages
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


#--------------
#Set the universal time
time = np.linspace(0, 1e-3,10000)


#-----------
#define the conditions

#define the values of the components in base units
Rs = 50         #Load resistence
R1 = 0.5        #Variable resistance circuit 2
R2 = 0.5        #Variable resistance circuit 1
R0 = 50         #Internal car resistance   
L1 = 27e-6      #Inductance of circuit 1
L2 = 27e-6      #Inductance of circuit 2
C1 = 1e-9       #Capacitance of circuit 1
C2 = 1e-9       #Capacitance of circuit 2

k = 0.3         #Coupling constant


#set the driving frequency
ω = 1e6


#set the voltage source, sinusoidal wave
imaginaryVs = np.exp(time*1j*ω)
realVs = np.real(imaginaryVs)


#define impdences from the components
#inductor impedence
L1imp = 1j*L1*ω
L2imp = 1j*L2*ω
Rsimp = Rs
R1imp = R1
R2imp = R2
R0imp = R0
C1imp = -1j/(ω*C1)
C2imp = -1j/(ω*C2)


    

#Variables calculated from the initial conditions
M = k*np.sqrt(L1*L2)                      #Mutual inductance
Z1 = Rsimp+R1imp+L1imp+C1imp              #Total inductance of circuit 1
Z2 = Rsimp+R1imp+L2imp+C2imp              #Total inductance of circuit 2



I1 = imaginaryVs/(Z1+(ω*M**2)/Z2)         #Current in circuit 1
I2 = 1j*ω*M*imaginaryVs/(Z1*Z2+(ω*M)**2)  #Current in circuit 2

powerIn  = imaginaryVs*I1
powerOut = I2**2*R0

Vout = Z2*I2

η = powerOut/powerIn


print("the efficiany is %.3f"%(np.max(η)))
print(np.format_float_scientific(np.max((np.real(I1))),precision=3))

     


#figure plotting
fig = plt.figure()
ax1 = plt.subplot(111)
#ax2 = plt.subplot(212)
ax1.grid()
plt.xlim(0,2e-4)
#ax2.grid()
ax1.plot(time,realVs)
ax1.plot(time,Vout)
plt.legend(['Vin',"Vout"])
#ax2.plot(time,I2,color="red")
