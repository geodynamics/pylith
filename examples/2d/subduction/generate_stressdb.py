#!/usr/bin/env python
"""
This script creates an initial stress spatial database from output stress
results.
"""

materials = ["concrust","oceancrust","conmantle","oceanmantle"]

# The code requires the numpy and h5py packages.
import numpy
import h5py

# Create coordinate system for spatial database
from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
from spatialdata.geocoords.CSCart import CSCart
cs = CSCart()
cs._configure()
cs.setSpaceDim(2)

for material in materials:

  filenameH5 = "output/step05-%s.h5" % material
  filenameDB = "grav_stress-%s.spatialdb" % material

  # Open HDF5 file and get coordinates, cells, and stress.
  h5 = h5py.File(filenameH5, "r")
  vertices = h5['geometry/vertices'][:]
  cells = numpy.array(h5['topology/cells'][:], dtype=numpy.int)
  stress = h5['cell_fields/stress'][0,:,:]
  if "mantleX" in material:
    stressZZ = h5['cell_fields/stress_zz_initial'][0,:,0]
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

  if "mantleX" in material:
    zeros = numpy.zeros(stressZZ.shape[0], dtype=numpy.float64)
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
               {'name': "viscous-strain-xy",
                'units': "None",
                'data': zeros},
               {'name': "viscous-strain-zz",
                'units': "None",
                'data': zeros},
          ]

  writer.write({'points': cellCenters,
                'coordsys': cs,
                'data_dim': 2,
                'values': values})

