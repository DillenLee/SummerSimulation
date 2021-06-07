import MutualInductance
import time
import numpy as np
import csv
import os
import datetime as dt


t0 = time.time()



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


name = 'TC_cy'
with open('/home/dillen/University/Python/Summer Project/Coil/%s.csv'%name,mode='w') as file:
    writer = csv.writer(file)
    writer.writerow([nt1, nt2, par1P, par1S, par2P, par2S, steps1, steps2, deltaK, deltaL])



#execute the loop


tstart = time.time()

for var in positionY:
    positionSecondary = [positionX,var,positionZ]
    Mi = MutualInductance.mutualInductance(typePrimary, typeSecondary, positionSecondary, nt1, nt2, par1P, par1S, par2P, par2S, steps1, steps2, deltaK, deltaL)
#    os.system('clear')

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


tEnd = time.time()

print(tEnd-t0)
