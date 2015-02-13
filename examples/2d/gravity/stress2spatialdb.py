#!/usr/bin/env python
"""
This script creates an initial stress spatial database from output stress
results.
"""

# The code requires the numpy and h5py packages.
import numpy
import h5py
import math
# import pdb
# pdb.set_trace()

# Define input/output files.
inFiles = ['output/grav_stress-elastic.h5',
           'output/grav_stress-viscoelastic.h5']
outFiles = ['grav_stress-elastic.spatialdb',
            'grav_stress-viscoelastic.spatialdb']
numFiles = len(inFiles)

# ----------------------------------------------------------------------
def convertStress(inFile, outFile):
  """
  Function to read HDF5 file and output stresses in a spatialdb.
  """

  # Open HDF5 file and get coordinates, cells, and stress.
  h5 = h5py.File(inFile, "r")
  coords = numpy.array(h5['geometry/vertices'][:], dtype=numpy.float64)
  cells = numpy.array(h5['topology/cells'][:], dtype=numpy.int)
  stress = numpy.array(h5['cell_fields/stress'][0,:,:], dtype=numpy.float64)
  numCells = cells.shape[0]
  h5.close()

  # Get cell centers for output.
  cellCoords = coords[cells,:]
  cellCenters = numpy.mean(cellCoords, axis=1)

  # Stack coordinates and stresses for output.
  outputArr = numpy.hstack((cellCenters, stress))

  # Write spatial database.
  dbHead1 = """// -*- C++ -*- (tell Emacs to use C++ mode)
//
// This spatial database gives the initial stresses for the model
//
#SPATIAL.ascii 1
SimpleDB {
  num-values = 3
  value-names =  stress-xx stress-yy stress-xy
  value-units =  Pa Pa Pa
  num-locs = %d
  """ % numCells

  dbHead2 = """data-dim = 2
  space-dim = 2
  cs-data = cartesian {
    to-meters = 1.0
    space-dim = 2
  }
}
"""
  f = open(outFile, 'w')
  f.write(dbHead1)
  f.write(dbHead2)
  numpy.savetxt(f, outputArr)
  f.close()

  return

# ======================================================================
# Loop over stress files.
for fileNum in range(numFiles):
  convertStress(inFiles[fileNum], outFiles[fileNum])
  
