#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 28 11:01:36 2021

@author: 
   ___  _ ____             __          
  / _ \(_) / /__ ___      / /  ___ ___ 
 / // / / / / -_) _ \    / /__/ -_) -_)
/____/_/_/_/\__/_//_/   /____/\__/\__/ 


Notes:
    Code adapted from 'Implementation of the Neumann Formula for
    Calculating the Mutual Inductance between Planar
    PCB Inductors' by C. L. W. Sonntag, E. A. Lomonova, and J. L. Duarte
    @TU/E
"""

import numpy as np
import matplotlib.pyplot as plt
import time 
import datetime as dt
import os
import csv




t0 = time.time()

#Some quick definitions
pi = np.pi

#Initial conditions
mu = 4*pi*1e-7              


#--------------
#Defined functions

#For part 1
#Single layered coil
def archimedeanSpiral(position,turns,innerRadius,outerRadius,steps):
    x0,y0,z0 = position
    lowerT = 0
    upperT = 2*pi*turns
    t = np.linspace(lowerT,upperT,steps)
    a = innerRadius
    b = (outerRadius-innerRadius)/(upperT)
    
    
    #The general form of the archimedian spiral in cylindrical coordinates is 
    #r = a + bθ, in this case it will be paremerised into its x and y coordinates
    #following x = r*cosθ and y = r*sinθ. The inner radius is detirmined by a
    #and the outer radius depends on the amount of turns according to the eqn
    #b = ΔD/(2π*t)
    
    x = x0+(np.cos(t)*(a+b*t))
    y = y0+(np.sin(t)*(a+b*t))
    z = z0
    xyzArray = np.array([x,y,z])
    
    pairedArray = []
    for i in range(len(t)):
        pairedArray.append([x[i],y[i],z])
    pairedArray = np.array(pairedArray)
    
    
    return xyzArray,pairedArray

#Multilayered coil
def coil(position,turns,radius,height,steps):
    x0,y0,z0 = position
    lowerT = 0
    upperT = 2*pi*turns
    t = np.linspace(lowerT,upperT,steps)
    r = radius
    
    #This is a simple coil converted from cylindrical coordinates to cartesian.
    
    x = x0+(r*np.cos(t))
    y = y0+(r*np.sin(t))
    z = z0+(height*(t/upperT))
    xyzArray = np.array([x,y,z])
    
    pairedArray = []
    for i in range(len(t)):
        pairedArray.append([x[i],y[i],z])
    pairedArray = np.array(pairedArray)

    
    return xyzArray,pairedArray

#For part 5
#creates the normalised vectors
def normalise(vector,ax,EuclidLength):
    n = []
    for i in range(1,ax):
        n.append((vector[i]-vector[i-1])/EuclidLength[i-1])
    return np.array(n)


#For part 4
#Calculates the magnitude of each vector pair
def magnitude(vector):
    a = []
    for i in range(1,len(vector)):
        difference = vector[i]-vector[i-1]
        mag = np.sqrt(np.dot(difference,difference))
        a.append(mag)
    return np.array(a)





#--------------

def mutualInductance(positionSecondary ,nt1 ,nt2 ,ir1 ,ir2 ,or1 ,or2 ,steps1 ,steps2 ,deltaK, deltaL):
    tstart = time.time()
    positionPrimary = [0,0,0]
    l = archimedeanSpiral(positionPrimary,nt1,ir1,or1,steps2)
    k = archimedeanSpiral(positionSecondary,nt2,ir2,or2,steps1)
    
    
    
    #--------------------
    #Part 1
    #Define the amount of line segments
    al = len(l[0][0])
    ak = len(k[0][0])
    
    #--------------------
    #Part 2
    #Define the amount of vertices
    #bl = al+1
    #bk = ak+1
    #NOT USED
    
    #--------------------
    #Part 3
    #The vector position of all the individual filaments
    ql = l[1]
    qk = k[1]
    
    #--------------------
    #Part 4
    #The length of each individual filament is found
    cl = magnitude(ql)
    ck = magnitude(qk)  
    
    nl = normalise(ql, al, cl)
    nk = normalise(qk, ak, ck)
    
    #Part 6
    
    
    constant = mu*deltaK*deltaL/(4*pi)
    M = 0
    for alpha in range(1,al):
        os.system('clear')
        print(alpha*100/al)
        tEnd = time.time()
        deltaT = tEnd-tstart
        tRemaining = deltaT*(1-alpha/al)/(1/al)
        tstart = tEnd
        print(dt.timedelta(seconds = tRemaining))
        for beta in range(1,ak):
            
            X = int(cl[alpha-1]/deltaL)
            E = int(ck[beta-1]/deltaK)
            
            for chi in range(0,X):
                for e in range(0,E):
                    #define some more variables
                    Kchi = ql[alpha-1] + chi*deltaL*nl[alpha-1]
                    Ke = qk[beta-1] + e*deltaK*nk[beta-1]
                    
                    
                    #do the actual calculation
                    #numerator, take the dot produc of the two
                    num = np.dot(nl[alpha-1],nk[beta-1])
                    #denominator, take the maginitude of the difference
                    diff = Kchi - Ke
                    den = np.sqrt(np.dot(diff,diff))
                    
                    #Sum up all the individual pieces
                    M += constant*num/den
                    
                    
                    
    
    fig = plt.figure()
    ax = plt.subplot(111,projection='3d')
    ax.plot(l[0][0],l[0][1],l[0][2])
    ax.plot(k[0][0],k[0][1],k[0][2])
                
    return M
    


    
#Do the loop thing
    
nt1 = nt2 = 20
ir1 = ir2 = 20 #mm
or1 = or2 = 40 #mm
steps1 = steps2 = 1000
deltaK = deltaL = 0.5
positionY = np.arange(-50,50,step = 5)  #mm

print(positionY)
'''
positionSecondary = [0,5,10]
mutualInductance(positionSecondary, nt1, nt2, ir1, ir2, or1, or2, steps1, steps2, deltaK, deltaL)
'''


mis = []
for var in positionY:
    positionSecondary = [0,var,10]
    Mi = mutualInductance(positionSecondary, nt1, nt2, ir1, ir2, or1, or2, steps1, steps2, deltaK, deltaL)   
    mis.append(Mi)
    with open('data.csv','a+',) as file:
        writer = csv.writer(file)
        writer.writerow([var,Mi])



print(time.time()-t0)
