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
import string




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
        pairedArray.append([x[i],y[i],z[i]])
    pairedArray = np.array(pairedArray)

    
    return xyzArray,pairedArray

#torus coil

def torus(position, turns ,ir,outr, steps):
    x0,y0,z0 = position
    t = np.linspace(0,steps,steps+1)
    R = ir
    r = outr - ir
    phi = 2*pi*t/steps
    theta = turns*phi
    
    x = x0+(R+r*np.cos(theta))*np.cos(phi)
    y = y0+(R+r*np.cos(theta))*np.sin(phi)
    z = z0+r*np.sin(theta)
    
    xyzArray = np.array([x,y,z])
    
    
    pairedArray = np.empty((len(x),3))

    for i in range(len(t)):
        pairedArray[i] =  [x[i],y[i],z[i]]
    
    return xyzArray, pairedArray


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

def mutualInductance(primaryType, secondaryType, positionSecondary ,ntP ,ntS ,Variable1P ,Variable1S ,Variable2P ,Variable2S ,steps1 ,steps2 ,deltaK, deltaL):
    positionPrimary = [0,0,0]
    if primaryType == 'S':
        l = archimedeanSpiral(positionPrimary,ntP,Variable1P,Variable2P,steps1)
    elif primaryType == 'C':
        l = coil(positionPrimary,ntP,Variable1P,Variable2P,steps1)
    elif primaryType == 'T':
        l = torus(positionPrimary,ntP,Variable1P,Variable2P,steps1)
    if secondaryType == 'S':     
        k = archimedeanSpiral(positionSecondary,ntS,Variable1S,Variable2S,steps2)
    elif secondaryType == 'C':
        k = coil(positionSecondary,ntS,Variable1S,Variable2S,steps2)
    elif secondaryType == 'T':
        k = torus(positionSecondary,ntS,Variable1S,Variable2S,steps2)
    
    
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
        
        print(alpha*100/al)
         
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
                    
                    
                    
    '''
    fig = plt.figure()
    ax = plt.subplot(111,projection='3d')
    ax.plot(l[0][0],l[0][1],l[0][2])
    ax.plot(k[0][0],k[0][1],k[0][2])
    ''' 
         
    return M
    



    
#Initial conditions

typePrimary = 'S'
typeSecondary = 'S' 
nt1 = 26
nt2 = 26
par1P = 5                   #mm
par1S = 5                   #mm  Inner radius for sprial and torus, radius for coil
par2P = 20                   #mm
par2S = 20                 #mm  Outer radius for spiral and torus, height for coil
steps1 = 100
steps2 = 100
deltaK = deltaL = .5
positionY = np.arange(-25,30,step = 5)  #mm
positionX = 0
positionZ = 100

# ---------------------
#To create a unique name for the file which contains all the info 
numbers = list(np.arange(0,65,1))
B64 = list(string.ascii_uppercase)+list(string.ascii_lowercase)+['0','1','2','3','4','5','6','7','8','9']+['+',"-"]

'''
def encode(message):
    aIndex = numbers.index(np.floor(message/64))
    bIndex = numbers.index(message%64)
    encode = B64[aIndex]+B64[bIndex]
    return encode

def decode(message):
    aIndex = B64.index(message[0])
    bIndex = B64.index(message[1])
    decode = numbers[aIndex]+numbers[bIndex]
    return decode

def numberOrletter(var,pos):
    if type(var) == int:
        return var
    else:
        return pos

turnP = encode(nt1)
turnS = encode(nt2)
var1P = encode(par1P)
var1S = encode(par1S)
var2P = encode(par2P)
var2S = encode(par2S)
middle = '-%s-%s-%s-'%(numberOrletter(positionX,'x'),numberOrletter(positionY,'y'),numberOrletter(positionZ,'z'))
name = typePrimary+typeSecondary+turnP+var1P+var2P+middle+turnS+var1S+var2S

print(name)
#---------------------

'''

name = 'TC'
with open('/home/dillen/University/Python/Summer Project/Coil/%s.csv'%name,mode='w') as file:
    writer = csv.writer(file)
    writer.writerow([nt1, nt2, par1P, par1S, par2P, par2S, steps1, steps2, deltaK, deltaL])
    


#execute the loop


tstart = time.time()

for var in positionY:
    positionSecondary = [positionX,var,positionZ]
    Mi = mutualInductance(typePrimary, typeSecondary, positionSecondary, nt1, nt2, par1P, par1S, par2P, par2S, steps1, steps2, deltaK, deltaL)   
    os.system('clear')
    varPosition = list(positionY).index(var)+1
    print(np.floor(varPosition*100/len(positionY)))
    tEnd = time.time()
    deltaT = tEnd-tstart
    tRemaining = deltaT*(1-varPosition/len(positionY))/(1/len(positionY))
    tstart = tEnd
    print(dt.timedelta(seconds = tRemaining))
       
    with open('/home/dillen/University/Python/Summer Project/Coil/%s.csv'%(name),'a') as file:
        writer = csv.writer(file)
        writer.writerow([var,Mi])



print(time.time()-t0)
