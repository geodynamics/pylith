import h5py

geom = "topo"

# Finite strain
sim = "%s_finite" % geom

filename = "output/%s-groundsurf.h5" % sim
h5 = h5py.File(filename, "r")
dispF = h5['vertex_fields/displacement'][:]
timeF = h5['time'][:,0,0]
h5.close()

filename = "output/%s-mantle.h5" % sim
h5 = h5py.File(filename, "r")
vstrainF = h5['cell_fields/viscous_strain'][:,0,:]
h5.close()


# Infinitesimal strain
sim = "%s_inf" % geom

filename = "output/%s-groundsurf.h5" % sim
h5 = h5py.File(filename, "r")
dispI = h5['vertex_fields/displacement'][:]
timeI = h5['time'][:,0,0]
h5.close()

filename = "output/%s-mantle.h5" % sim
h5 = h5py.File(filename, "r")
vstrainI = h5['cell_fields/viscous_strain'][:,0,:]
h5.close()

from pyre.units.time import year
timeF /= year.value
timeI /= year.value

import pylab
iPt = 1
pylab.subplot(2,1,1)
pylab.plot(timeF, dispF[:,iPt,1], 'r-', timeI, dispI[:,iPt,1], 'b--')

pylab.subplot(2,1,2)
pylab.plot(timeF, vstrainF[:,:], 'r-', timeI, vstrainI[:,:], 'b--')

pylab.show()

