#!/usr/bin/env python
"""
This script creates a spatial database for the initial stress and state
variables for a Maxwell plane strain material.
"""

sim = "topo-finite"

import numpy
import h5py

from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
from spatialdata.geocoords.CSCart import CSCart
cs = CSCart()
cs._configure()
cs.setSpaceDim(2)

filenameH5 = "output/%s-statevars.h5" % sim
filenameDB = "%s-restart-initialstate.spatialdb" % sim

# Open HDF5 file and get coordinates, cells, and stress.
h5 = h5py.File(filenameH5, "r")
vertices = h5['geometry/vertices'][:]
cells = numpy.array(h5['topology/cells'][:], dtype=numpy.int)
stress = h5['cell_fields/stress'][0,:,:]
h5.close()

# Get cell centers for output.
cellCoords = vertices[cells,:]
qpts = numpy.array([[ 0.62200847,  0.16666667,  0.0446582,   0.16666667],
                    [ 0.16666667,  0.62200847,  0.16666667,  0.0446582 ],
                    [ 0.16666667,  0.0446582,   0.16666667,  0.62200847],
                    [ 0.0446582,   0.16666667,  0.62200847,  0.16666667]], dtype=numpy.float64)
print cellCoords.shape
qptsX = numpy.dot(cellCoords[:,:,0], qpts)
qptsY = numpy.dot(cellCoords[:,:,1], qpts)

print qptsX.shape
cellXY = numpy.hstack((qptsX, qptsY))
print cellXY.shape
cellCenters = numpy.mean(cellCoords, axis=1)

# Create writer for spatial database file
writer = SimpleIOAscii()
writer.inventory.filename = filenameDB
writer._configure()

values = [{'name': "stress-xx",
           'units': "Pa",
           'data': stress[:,0]},
          {'name': "stress-yy",
           'units': "Pa",
           'data': stress[:,1]},
          {'name': "stress-xy",
           'units': "Pa",
           'data': stress[:,2]},
        ]

stressZZ = stress[:,0]
zeros = numpy.zeros(stressZZ.shape)
values += [{'name': "stress-zz-initial",
            'units': "Pa",
            'data': stressZZ},
           {'name': "total-strain-xx",
            'units': "None",
            'data': zeros},
           {'name': "total-strain-yy",
            'units': "None",
            'data': zeros},
           {'name': "total-strain-xy",
            'units': "None",
            'data': zeros},
           
           {'name': "viscous-strain-xx",
            'units': "None",
            'data': zeros},
           {'name': "viscous-strain-yy",
            'units': "None",
            'data': zeros},
           {'name': "viscous-strain-zz",
            'units': "None",
            'data': zeros},
           {'name': "viscous-strain-xy",
            'units': "None",
            'data': zeros},
         ]

writer.write({'points': cellCenters,
              'coordsys': cs,
              'data_dim': 1,
              'values': values})

# End of file
