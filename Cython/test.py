import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import box

t0 = time.time()
data = box.torus([0,0,0],260,4.3,8.2,5000)[0]
print(time.time()-t0)

fig = plt.figure()

shape = [-10,10]


ax = Axes3D(fig,xlim=shape,ylim=shape,zlim=shape)
ax.plot(data[0],data[1],data[2])
plt.savefig('torus.png',dpi = 600)
