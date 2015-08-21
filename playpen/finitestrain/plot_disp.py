import h5py

sim = "topo"


filename = "output/%s-finite-domain.h5" % sim
h5 = h5py.File(filename, "r")
dispF = h5['vertex_fields/displacement'][:]
timeF = h5['time'][:,0,0]
h5.close()

if True:
    filename = "output/%s-inf-domain.h5" % sim
    h5 = h5py.File(filename, "r")
    dispI = h5['vertex_fields/displacement'][:]
    timeI = h5['time'][:,0,0]
    h5.close()


print "Finite",dispF[-1,:,:]
print "Infinitesimal",dispI[-1,:,:]


import pylab
iPt = 1
pylab.plot(timeF, dispF[:,iPt,1], 'r-', timeI, dispI[:,iPt,1], 'b--')
pylab.show()

