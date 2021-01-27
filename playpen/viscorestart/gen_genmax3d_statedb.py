#!/usr/bin/env python
"""
This script creates a spatial database for the initial stress and state
variables for a Generalized Maxwell 3D material.
"""

material = "genmax3d"

import numpy
import h5py

from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
from spatialdata.geocoords.CSCart import CSCart
cs = CSCart()
cs._configure()
cs.setSpaceDim(3)

filenameH5 = "output/grav_static_%s-visco.h5" % material
filenameDB = "grav_statevars-%s.spatialdb" % material

# Open HDF5 file and get coordinates, cells, and stress.
h5 = h5py.File(filenameH5, "r")
vertices = h5['geometry/vertices'][:]
cells = numpy.array(h5['topology/cells'][:], dtype=numpy.int)
stress = h5['cell_fields/stress'][0,:,:]
strain = h5['cell_fields/total_strain'][0,:,:]
strainViscous1 = h5['cell_fields/viscous_strain_1'][0,:,:]
strainViscous2 = h5['cell_fields/viscous_strain_2'][0,:,:]
strainViscous3 = h5['cell_fields/viscous_strain_3'][0,:,:]
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
           'data': stress[:, 0]},
          {'name': "stress-yy",
           'units': "Pa",
           'data': stress[:, 1]},
          {'name': "stress-zz",
           'units': "Pa",
           'data': stress[:, 2]},
          {'name': "stress-xy",
           'units': "Pa",
           'data': stress[:, 3]},
          {'name': "stress-yz",
           'units': "Pa",
           'data': stress[:, 4]},
          {'name': "stress-xz",
           'units': "Pa",
           'data': stress[:, 5]},
        ]

#if "mantle" in material:
if True:
  zeros = numpy.zeros(stress[:, 0].shape)
  values += [{'name': "total-strain-xx",
              'units': "None",
              'data': zeros},
             {'name': "total-strain-yy",
              'units': "None",
              'data': zeros},
             {'name': "total-strain-zz",
              'units': "None",
              'data': zeros},
             {'name': "total-strain-xy",
              'units': "None",
              'data': zeros},
             {'name': "total-strain-yz",
              'units': "None",
              'data': zeros},
             {'name': "total-strain-xz",
              'units': "None",
              'data': zeros},

             {'name': "viscous-strain-1-xx",
              'units': "None",
              'data': strainViscous1[:, 0]},
             {'name': "viscous-strain-1-yy",
              'units': "None",
              'data': strainViscous1[:, 1]},
             {'name': "viscous-strain-1-zz",
              'units': "None",
              'data': strainViscous1[:, 2]},
             {'name': "viscous-strain-1-xy",
              'units': "None",
              'data': strainViscous1[:, 3]},
             {'name': "viscous-strain-1-yz",
              'units': "None",
              'data': strainViscous1[:, 4]},
             {'name': "viscous-strain-1-xz",
              'units': "None",
              'data': strainViscous1[:, 5]},

             {'name': "viscous-strain-2-xx",
              'units': "None",
              'data': strainViscous2[:, 0]},
             {'name': "viscous-strain-2-yy",
              'units': "None",
              'data': strainViscous2[:, 1]},
             {'name': "viscous-strain-2-zz",
              'units': "None",
              'data': strainViscous2[:, 2]},
             {'name': "viscous-strain-2-xy",
              'units': "None",
              'data': strainViscous2[:, 3]},
             {'name': "viscous-strain-2-yz",
              'units': "None",
              'data': strainViscous2[:, 4]},
             {'name': "viscous-strain-2-xz",
              'units': "None",
              'data': strainViscous2[:, 5]},

             {'name': "viscous-strain-3-xx",
              'units': "None",
              'data': strainViscous3[:, 0]},
             {'name': "viscous-strain-3-yy",
              'units': "None",
              'data': strainViscous3[:, 1]},
             {'name': "viscous-strain-3-zz",
              'units': "None",
              'data': strainViscous3[:, 2]},
             {'name': "viscous-strain-3-xy",
              'units': "None",
              'data': strainViscous3[:, 3]},
             {'name': "viscous-strain-3-yz",
              'units': "None",
              'data': strainViscous3[:, 4]},
             {'name': "viscous-strain-3-xz",
              'units': "None",
              'data': strainViscous3[:, 5]},
        ]

writer.write({'points': cellCenters,
              'coordsys': cs,
              'data_dim': 0,
              'values': values})

# End of file
