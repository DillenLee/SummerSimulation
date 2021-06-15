#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author:
   ___  _ ____             __
  / _ \(_) / /__ ___      / /  ___ ___
 / // / / / / -_) _ \    / /__/ -_) -_)
/____/_/_/_/\__/_//_/   /____/\__/\__/

"""
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from numpy import pi
from mpl_toolkits.mplot3d import Axes3D
import time
import csv

t0 = time.time()

# The first part of this code is to implement a discrete form of Biot-Savart law

# Define constants
mu = 4*pi*1e-7

# Initial conditions
I = 1/300           # (A) current through the transmitting coil
radiusTrans = 47.5  # (mm) radius of transmitting coil
radiusRec   = 23.5  # (mm) radius of receiving coil
Nr = 18           # amount of turns in receiving coil
Nt = 10.5           # amount of turns in transmitting coil
zPos =  23          # (mm) Z distance between coils
tHeight = 4.75       # (mm) The height of the transmitting coil

# The simplest case would be two circles so this function will create an array
# discretising the circle
def circle(position,r,steps):
    x0,y0,z0 = position
    theta = np.linspace(0,2*pi,steps)
    x = x0+r*np.cos(theta)
    y = y0+r*np.sin(theta)
    z = np.full(len(theta),z0)
    array = np.array([x,y,z])
    return array

# Needed later

def archimedeanSpiral(position,turns,radius,steps):
    x0,y0,z0 = position
    lowerT = 0
    upperT = 2*pi*turns
    t = np.linspace(lowerT,upperT,steps)
    b = radius/(upperT)

    #The general form of the archimedian spiral in cylindrical coordinates is
    #r = a + bθ, in this case it will be paremerised into its x and y coordinates
    #following x = r*cosθ and y = r*sinθ. The inner radius is detirmined by a
    #and the outer radius depends on the amount of turns according to the eqn
    #b = ΔD/(2π*t)

    x = x0+(np.cos(t)*b*t)
    y = y0+(np.sin(t)*b*t)
    z = z0

    pairedArray = np.ndarray((len(t),3))
    for i in range(len(t)):
        pairedArray[i] = [x[i],y[i],z]

    return pairedArray
    # return np.array([x,y,z])

#Multilayered coil, Needed later too
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

    return xyzArray






# This function simply returns a normalised a vector
def hat(vector):
    magnitude = np.sqrt(vector.dot(vector))
    vectorHat = vector/magnitude
    return vectorHat

# Now here comes the actual meat of the code, the Biot-Savart Law
# Object has array type with each dimension as an individual array i.e. all the xs in one array within 'object'
# Position is also an array of the x,y,z (r) position of the point to be calculated
def BiotSavart(object,position):
    x,y,z = object
    dl = np.empty((len(x),3))           # Make an array which contains the direction vector
    for i in range(1,len(x)):           # of the position vectors in the object (dl)
        xt = x[i]-x[i-1]
        yt = y[i]-y[i-1]
        zt = z[i]-z[i-1]
        dl[i] = np.array([xt,yt,zt])

    B = np.ndarray(3)
    for i in range(len(dl)):                                # Now we do the actual sum
        rPrime = position-np.array([x[i],y[i],z[i]])        # rPrime is the distance between the position and the coil position element dl
        rPrimeHat = hat(rPrime)
        cross = np.cross(dl[i],rPrimeHat)
        toAdd = cross/rPrime.dot(rPrime)
        B += toAdd
        # print(B)
        # time.sleep(0)

    B *= (mu*I)/(4*pi)

    return B


# Create a file for the Data
name = 'BSsmallR'
with open('/home/dillen/University/Python/SummerSimulation/Data/%s.csv'%name,mode='w') as file:
    writer = csv.writer(file)
    writer.writerow(['Distance (mm)','Mutual Inductance (mH)','Small coil'])



# We loop over a bunch of y distances to get the mutual inductance for different values



for y in range(-225,235,5):
    # Now that we can give a B vector to every point in space time to find the magnetic flux
    # We can model the area of the solenoid as a simple circle.
    # First let's determine how many magnetic field points we want and create a grid
    # which will be the positions. An archimedian spiral from the centre outwards will
    # be a pretty good approximation to cover a circle shape
    rPositions = archimedeanSpiral([0,y,zPos+tHeight],200,radiusRec,1000)

    # Now calculate the B field at every position and take an average value
    # Also since we know that only the z component will matter regarding flux
    # lets just only take that position for now. (The dS is facing purely upwards)
    Btotal = 0
    for pos in rPositions:
        # Bi = BiotSavart(circle([0,0,23],radiusRec,1000),pos)
        Bi = BiotSavart(coil([0,0,0],Nt,radiusTrans,tHeight,1000),pos)
        Btotal += Bi[2]                                 #Take the z component and sum
    # The average B is the total B divided by the sampling number
    Bavg = Btotal/len(rPositions)

    #flux is simply BA (no cos needed as already perpendicular)
    phi = Bavg*pi*(radiusRec)**2

    #now we can find the magnetic inductance using the simple formula M = N*phi/I
    M = Nr*phi/I

    with open('/home/dillen/University/Python/SummerSimulation/Data/%s.csv'%(name),'a') as file:
        writer = csv.writer(file)
        writer.writerow([y,M])
    print(y*100/500)

print(time.time()-t0)
