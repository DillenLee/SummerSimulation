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

import matplotlib.pyplot as plt
import numpy as np


def extract(ID,rowsToSkip):
    x,y = np.loadtxt(ID, skiprows=rowsToSkip, delimiter=',', unpack=True,usecols=[0,1])
    return x,y

BiotSavartD, BiotSavartM = extract('Data/BSlargeR.csv', 1)
compD,compM = extract('Data/R1.csv', 1)
BSDsmall,BSMsmall = extract('Data/BSsmallR.csv',1)
compDsmall,compMsmall = extract('Data/R2.csv',1)

fig = plt.figure()
ax = plt.subplot(111)
ax.grid()
ax.set_xlabel('Distance from equilibrium (mm)')
ax.set_ylabel('Mutual inductance (mH)')
ax.plot(BiotSavartD,BiotSavartM,label = 'Biot Savart method large')
ax.plot(compD,compM,label = 'Neumann method large')
ax.plot(BSDsmall,BSMsmall,label='Biot Savart method large')
ax.plot(compDsmall,compMsmall,label='Neumann method small')
ax.legend(loc = 'upper right')

# plt.show()
plt.savefig('NeumannVsBiotSavart.png',dpi=300)
