#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

## @file unittests/libtests/feassemble/data/Mesh1DLinear.py

## @brief Python container for data associated with finite-element
## mesh in 1-D space with linear basis funtions.

from IntegratorInertia import IntegratorInertia

import numpy

# ----------------------------------------------------------------------

# Mesh1DLinear class
class Mesh1DLinear(object):
  """
  Python application for generating C++ data files for testing C++
  IntegratorInertia object with 1-D cell and linear basis functions.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self):
    """
    Constructor.
    """
    # Mesh information
    self.spaceDim = 1
    self.cellDim = 1
    self.numVertices = 2
    self.numCells = 1
    self.vertices = numpy.array( [[-0.25], [2.0]], dtype=numpy.float64)
    self.cells = numpy.array( [[0, 1]], dtype=numpy.int32)

    # Quadrature information
    from Quadrature1DLinear import Quadrature1DLinear
    self.quadrature = Quadrature1DLinear()

    # Material information
    self.matType = "ElasticStress1D"
    self.matDBFilename = "data/elasticstress1d.spatialdb"
    self.matId = 0
    self.label = "elastic 1-D"
    self.mu = 3.0e+10
    self.lambda = self.mu
    self.density = 2500.0

    # Input fields
    self.dt = 0.01
    self.fieldTpdt = numpy.array( [[1.2], [1.7]], dtype=numpy.float64)
    self.fieldT = numpy.array( [[1.1], [1.5]], dtype=numpy.float64)
    self.fieldTmdt = numpy.array( [[1.0], [1.3]], dtype=numpy.float64)
    return
  

# End of file 
