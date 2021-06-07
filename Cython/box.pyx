import numpy as np
cimport numpy as np

pi = np.pi

cpdef tuple torus(list position, double turns ,double ir, double outr,int steps):
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
