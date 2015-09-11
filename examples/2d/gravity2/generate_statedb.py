#!/usr/bin/env python
"""
This script creates a spatial database for the initial stress and state
variables for a Maxwell plane strain material.
"""

materials = ["crust","mantle"]

import numpy
import h5py

from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
from spatialdata.geocoords.CSCart import CSCart
cs = CSCart()
cs._configure()
cs.setSpaceDim(2)

# Basis functions for quad4 cell evaluated at quadrature points. Use
# to compute coordinate of quadrature points in each cell from
# coordinates of vertices.
qpts = numpy.array([[ 0.62200847,  0.16666667,  0.0446582,   0.16666667],
                    [ 0.16666667,  0.62200847,  0.16666667,  0.0446582 ],
                    [ 0.16666667,  0.0446582,   0.16666667,  0.62200847],
                    [ 0.0446582,   0.16666667,  0.62200847,  0.16666667]], dtype=numpy.float64)


def calcQuadCoords(vertices, cells, qpts):
  """Compute coordinates of quadrature points."""
  nqpts = qpts.shape[0]
  ncells = cells.shape[0]
  spaceDim = vertices.shape[1]

  quadCoords = numpy.zeros((ncells, nqpts, spaceDim), dtype=numpy.float64)
  cellCoords = vertices[cells,:]
  for iDim in xrange(spaceDim):
    quadCoords[:,:,iDim] = numpy.dot(cellCoords[:,:,iDim], qpts)

  quadCoords = quadCoords.reshape((ncells*nqpts, spaceDim))
  return quadCoords


for material in materials:

  filenameH5 = "output/topo_inf-%s.h5" % material
  filenameDB = "grav_statevars-%s.spatialdb" % material

  # Open HDF5 file and get coordinates, cells, and stress.
  h5 = h5py.File(filenameH5, "r")
  vertices = h5['geometry/vertices'][:]
  cells = numpy.array(h5['topology/cells'][:], dtype=numpy.int)
  stress = h5['cell_fields/cauchy_stress'][-1,:,:]
  if "mantle" in material:
    vstrain = h5['cell_fields/viscous_strain'][-1,:,:]
  h5.close()
  
  # Compute coordinates of quadrature points.
  quadCoords = calcQuadCoords(vertices, cells, qpts)
  
  nqpts = qpts.shape[0]
  ncells = cells.shape[0]
  nvalues = stress.shape[1]/nqpts
  stress = stress.reshape((ncells*nqpts, nvalues))

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
  
  if "mantle" in material:
    nvalues = vstrain.shape[1]/nqpts
    vstrain = vstrain.reshape((ncells*nqpts, nvalues))
    stressZZ = 0.5*(stress[:,0] + stress[:,1])
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
                'data': vstrain[:,0]},
               {'name': "viscous-strain-yy",
                'units': "None",
                'data': vstrain[:,1]},
               {'name': "viscous-strain-zz",
                'units': "None",
                'data': vstrain[:,2]},
               {'name': "viscous-strain-xy",
                'units': "None",
                'data': vstrain[:,3]},
             ]

  writer.write({'points': quadCoords,
                'coordsys': cs,
                'data_dim': 2,
                'values': values})

# End of file
