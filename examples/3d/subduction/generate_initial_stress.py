#!/usr/bin/env nemesis
"""
This script creates a spatial database for the initial stress using results
from a previous PyLith run.
"""

dbPrefix = "initial_stress"
matPrefix = "step08a"
materials = ["crust","mantle","slab","wedge"]

import numpy
import h5py
# import pdb
# pdb.set_trace()

from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
from mesh.coordsys import cs_mesh
cs = cs_mesh()

def getCellCenters(vertices, cells):
  """
  Function to compute cell centers.
  """
  cellCoords = vertices[cells, :]
  cellCenters = numpy.mean(cellCoords, axis=1)

  return cellCenters


for material in materials:

  filenameH5 = "output/%s-%s.h5" % (matPrefix, material)
  filenameDB = "spatialdb/%s-%s.spatialdb" % (dbPrefix, material)

  # Open HDF5 file and get coordinates, cells, and stress.
  h5 = h5py.File(filenameH5, "r")
  vertices = h5['geometry/vertices'][:]
  cells = numpy.array(h5['topology/cells'][:], dtype=numpy.int)
  stress = h5['cell_fields/stress'][0,:,:]
  h5.close()
  
  # Compute coordinates of quadrature points.
  quadCoords = getCellCenters(vertices, cells)
  
  ncells = cells.shape[0]
  nvalues = stress.shape[1]

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
            {'name': "stress-zz",
             'units': "Pa",
             'data': stress[:,2]},
            {'name': "stress-xy",
             'units': "Pa",
             'data': stress[:,3]},
            {'name': "stress-yz",
             'units': "Pa",
             'data': stress[:,4]},
            {'name': "stress-xz",
             'units': "Pa",
             'data': stress[:,5]},
          ]
  
  writer.write({'points': quadCoords,
                'coordsys': cs,
                'data_dim': 3,
                'values': values})

# End of file
