#!/usr/bin/env python
"""
This script creates a spatial database for the initial stress and state
variables for a Maxwell plane strain material.
"""

material = "maxps-oceanmantle"

import numpy
import h5py

from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
from spatialdata.geocoords.CSCart import CSCart
cs = CSCart()
cs._configure()
cs.setSpaceDim(2)

filenameH5 = "output/grav_static_%s.h5" % material
filenameDB = "grav_statevars-%s.spatialdb" % material

# Open HDF5 file and get coordinates, cells, and stress.
h5 = h5py.File(filenameH5, "r")
vertices = h5['geometry/vertices'][:]
cells = numpy.array(h5['topology/cells'][:], dtype=numpy.int)
stress = h5['cell_fields/stress'][0,:,:]
strain = h5['cell_fields/total_strain'][0,:,:]
strainViscous = h5['cell_fields/viscous_strain'][0,:,:]
h5.close()

# Get cell centers for output.
cellCoords = vertices[cells,:]
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

#if "mantle" in material:
if True:
  density = 4000.0
  vs = 5600.0
  vp = 10000.0
  modulus_mu = density*vs**2
  modulus_lambda = density*vp**2 - 2.0*modulus_mu
  stressZZ = modulus_lambda * (strain[:,0] + strain[:,1])
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
              'data': strainViscous[:,0]},
             {'name': "viscous-strain-yy",
              'units': "None",
              'data': strainViscous[:,1]},
             {'name': "viscous-strain-zz",
              'units': "None",
              'data': strainViscous[:,2]},
             {'name': "viscous-strain-xy",
              'units': "None",
              'data': strainViscous[:,3]},
        ]

writer.write({'points': cellCenters,
              'coordsys': cs,
              'data_dim': 0,
              'values': values})

# End of file
