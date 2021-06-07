#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  5 19:31:49 2021

@author:
   ___  _ ____             __
  / _ \(_) / /__ ___      / /  ___ ___
 / // / / / / -_) _ \    / /__/ -_) -_)
/____/_/_/_/\__/_//_/   /____/\__/\__/

"""

#import the necessary packages
import numpy as np
import matplotlib.pyplot as plt

#define a function to easily extract the x and y components
def extract(ID,rowsToSkip):
    x,y = np.loadtxt(ID, skiprows=rowsToSkip, delimiter=',', unpack=True,usecols=[0,1])
    return x,y

#--------------------------
#Compare the large coil radius 38.5mm

#extract the computational data points
compD,compM = extract('Data/R1.csv', 2)     #compD is distance (mm), compM is mutual inductance (mH)

#extract the experimental data
expT, expE = extract('Data/exp1.csv',3)     #expT is time (s), expE is induced EMF ε (mV)

#This part will compare the mutual inductance between the experimental
#and the computational simulation

#---------------------------

#Initial conditions
Vs = 2              # (V) Source, driving potential difference
Rs = 0.5            # (Ω) Source resistance
f = 10814.917055/(2*np.pi)          # (Hz) Source, driving frequency
velocity = 0.05     # (m/s) lift velocity




#----calculations-----
Is = Vs/Rs          # (A) Driving current
ω = 2*np.pi*f       # (rad/s) Source driving angular frequency




#---------------------------
#Convert the time scale of the experimental to distance scale,
#assume constant velocity, Also centre at x = 0 with largest EMF value


expD = velocity*expT*1e3            #since the computational distance is in mm convert experimental from m to mm

maxE = np.amax(expE)                #take the largest EMF
maxEIndex = np.where(expE == maxE)  #and find the index position
expD -= expD[maxEIndex]             #and now subtract by that distance


expM = expE/(Is*ω*np.cos(ω*expD/(velocity*1e3)))                     #Mutual inductance equation for t = n*ω/2π
                                    #not comletely correct as there should be a
                                    #sin(ωt) component in the denominator but it
                                    #kinda fucks up with the discretisation of the
                                    #picoscope

#clean the broken data points
deletePoints = []
for i in range(len(expM)):
    if np.abs(expM[i]) > np.amax(compM)*20:
        deletePoints.append(i)

expM = np.delete(expM,deletePoints)
expD = np.delete(expD,deletePoints)


#----plot the data----

#Control parameters for inspection
minX = -200
maxX = 200


#boring matplotlib stuff
fig = plt.figure(dpi=200)
ax = plt.subplot(111)
ax.set_title('Comparison of mutual induction')
ax.grid()
ax.set_xlim([minX,maxX])
ax.set_xlabel('Distance from equilibrium (mm)')
ax.set_ylabel('Mutual inductance (mH)')
ax.plot(expD,expM)
ax.plot(compD,compM)
ax.legend(['Experimental data','Simulation model'],loc ='lower right')
fig.savefig('MutualInductance.png',dpi=600,bbox_inches='tight')


#----------------------------------------------
#Comparison of the coupling coefficient K for the large coil

#define the coupling equation
def coupling(M,L1,L2):
    return M/np.sqrt(L1*L2)

#set some variables in case it changes but I doubt it
L1 = 26e-6                  # (H) inductance of primary coil
L2 = 26e-6                  # (H) inductance of the secondary coil

#Define the coupling coefficents
compK = coupling(compM*1e-3,L1,L2)       #for the computational model with conversion from mH to H
expK = coupling(expM*1e-3,L1,L2)         #for the experimental model with conversion from mH to H

#Control parameters for inspection
minX = -200
maxX = 200

fig2 = plt.figure()
ax2 = plt.subplot(111)
ax2.set_title('Comparison of coupling coefficent')
ax2.grid()
ax2.set_xlim([minX,maxX])
ax2.set_xlabel('Distance from equilibrium (mm)')
ax2.set_ylabel('Coupling coefficient, K')
ax2.plot(expD,expK,color='green')
ax2.plot(compD,compK,color='red')
ax2.legend(['Experimental data','Simulation model'],loc ='upper right')
fig2.savefig('CouplingCoefficient.png',dpi=600,bbox_inches='tight')
plt.tight_layout()
plt.show()
