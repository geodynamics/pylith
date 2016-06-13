#!/usr/bin/env python
"""
This script creates a spatial database for the initial stress and state
variables for a Maxwell plane strain material.
"""

sim = "gravity_vardensity"
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
# coordinates of vertices. Note the order must correspond to the order
# of the data at the quadrature points in the output.
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
    quadCoords[:,:,iDim] = numpy.dot(cellCoords[:,:,iDim], qpts.transpose())

  quadCoords = quadCoords.reshape((ncells*nqpts, spaceDim))
  return quadCoords

for material in materials:

  filenameH5 = "output/%s-%s.h5" % (sim, material)
  filenameDB = "%s_statevars-%s.spatialdb" % (sim, material)

  # Open HDF5 file and get coordinates, cells, and stress.
  h5 = h5py.File(filenameH5, "r")
  vertices = h5['geometry/vertices'][:]
  tindex = -1
  cells = numpy.array(h5['topology/cells'][:], dtype=numpy.int)
  stress = h5['cell_fields/stress'][tindex,:,:]
  if "mantle" in material:
    vstrain = h5['cell_fields/viscous_strain'][tindex,:,:]
  h5.close()
  
  # Compute coordinates of quadrature points.
  quadCoords = calcQuadCoords(vertices, cells, qpts)
  
  nqpts = qpts.shape[0]
  ncells = cells.shape[0]
  nvalues = stress.shape[1]/nqpts  

  # Check to make sure output included all quadrature points (CellFilterAvg was not used).
  if stress.shape[1] == 3:
    raise ValueError("Found %d stress values for each cell. Expected 12 stress values (stress_xx, stress_yy, and stress_xy at 4 quadrature points) for each cell. Turn off CellFilterAvg in pylithapp.cfg." % stress.shape[1])

  if stress.shape[1] != nqpts*3:
    raise ValueError("Found %d stress values for each cell. Expected 12 stress values (stress_xx, stress_yy, and stress_xy at 4 quadrature points) for each cell. Did you turn off CellFilterAvg in pylithapp.cfg?" % stress.shape[1])

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
    stressZZ = 0.5*(stress[:,0]+stress[:,1])
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
