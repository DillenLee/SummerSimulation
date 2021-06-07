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
cimport numpy as np
cimport cython




#Some quick definitions
cdef double pi = np.pi

#Initial conditions
cdef double mu = 4*pi*1e-7


#--------------
#Defined functions

#For part 1
#Single layered coil
@cython.boundscheck(False)
cpdef inline tuple archimedeanSpiral(list position,double turns,double innerRadius,double outerRadius,int steps):
    cdef double x0 = position[0]
    cdef double y0 = position[1]
    cdef double z0 = position[2]
    cdef np.ndarray x
    cdef np.ndarray y
    cdef double z
    cdef np.ndarray t
    cdef np.ndarray xyzArray
    cdef int i
    cdef double a
    cdef double b
    cdef upperT = 2*pi*turns

    t = np.linspace(0,upperT,steps)
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

    pairedArray = np.empty((t.size,3))

    for i in np.arange(t.size):
        pairedArray[i] =  [x[i],y[i],z]

    return xyzArray,pairedArray

#Multilayered coil
@cython.boundscheck(False)
cpdef inline tuple coil(list position,double turns, double radius,double height,int steps):
    cdef double x0
    cdef double y0
    cdef double z0
    cdef np.ndarray x
    cdef np.ndarray y
    cdef np.ndarray z
    cdef np.ndarray t
    cdef np.ndarray xyzArray
    cdef int i
    cdef double r



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

    pairedArray = np.empty((t.size,3))

    for i in np.arange(t.size):
        pairedArray[i] =  [x[i],y[i],z[i]]


    return xyzArray,pairedArray

#torus coil
@cython.boundscheck(False)
cpdef inline tuple torus(list position, double turns ,double ir, double outr,int steps):
    cdef double x0
    cdef double y0
    cdef double z0
    cdef np.ndarray t
    cdef double R
    cdef np.ndarray phi
    cdef double r
    cdef np.ndarray theta
    cdef np.ndarray x
    cdef np.ndarray y
    cdef np.ndarray z
    cdef int i
    cdef np.ndarray pairedArray

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


    pairedArray = np.empty((t.size,3))

    for i in np.arange(t.size):
        pairedArray[i] =  [x[i],y[i],z[i]]

    return xyzArray,pairedArray


#For part 5
#creates the normalised vectors
@cython.boundscheck(False)
cpdef inline np.ndarray normalise(np.ndarray vector,int ax,np.ndarray EuclidLength):
    cdef np.ndarray n = np.empty((ax,3))
    cdef int i

    for i in np.arange(1,ax):
        n[i-1] = ((vector[i]-vector[i-1])/EuclidLength[i-1])
    return n


#For part 4
#Calculates the magnitude of each vector pair
@cython.boundscheck(False)
cpdef inline np.ndarray magnitude(np.ndarray vector):
    cdef np.ndarray a = np.empty((vector.size))
    cdef np.ndarray difference
    cdef int i
    for i in np.arange(1,len(vector)):
        difference = vector[i]-vector[i-1]
        a[i-1] = np.sqrt(difference.dot(difference))
    return a






#--------------

@cython.boundscheck(False)
#cpdef double mutualInductance(str primaryType='T',str secondaryType='T',list positionSecondary=[10,-25,100] ,double ntP=26.0  ,double ntS=26.0 ,double Variable1P=5.0 ,double Variable1S=5.0 ,double Variable2P=20.0 ,double Variable2S=20.0 ,int steps1=100 ,int steps2=100 ,double deltaK=.5,double deltaL=.5):
cpdef double mutualInductance(str primaryType,str secondaryType,list positionSecondary ,double ntP ,double ntS ,double Variable1P ,double Variable1S ,double Variable2P ,double Variable2S ,int steps1 ,int steps2 ,double deltaK,double deltaL):
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
    cdef int al
    cdef int ak
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
    cdef np.ndarray ql
    cdef np.ndarray qk
    ql = l[1]
    qk = k[1]

    #--------------------
    #Part 4
    #The length of each individual filament is found
    cdef np.ndarray cl
    cdef np.ndarray ck

    cl = magnitude(ql)
    ck = magnitude(qk)

    cdef np.ndarray nl
    cdef np.ndarray nk

    nl = normalise(ql, al, cl)
    nk = normalise(qk, ak, ck)

    #Part 6

    cdef double constant
    cdef double M
    cdef int alpha
    cdef int beta
    cdef int chi
    cdef int e
    cdef int X
    cdef int E
    cdef np.ndarray Kchi
    cdef np.ndarray Ke
    cdef np.ndarray diff
    cdef double num
    cdef double den




    constant = mu*deltaK*deltaL/(4*pi)

    M = 0.0
    for alpha in np.arange(1,al):

        print(alpha*100/al)

        for beta in np.arange(1,ak):

            X = cl[alpha-1]/deltaL
            E = ck[beta-1]/deltaK

            for chi in np.arange(X):
                for e in np.arange(E):
                    #define some more variables
                    Kchi = ql[alpha-1] + chi*deltaL*nl[alpha-1]
                    Ke = qk[beta-1] + e*deltaK*nk[beta-1]

                    #do the actual calculation
                    #numerator, take the dot produc of the two
                    num = nl[alpha-1].dot(nk[beta-1])
                    #denominator, take the maginitude of the difference
                    diff = Kchi - Ke
                    den = np.sqrt(diff.dot(diff))


                    #Sum up all the individual pieces
                    M += constant*num/den

    return M
