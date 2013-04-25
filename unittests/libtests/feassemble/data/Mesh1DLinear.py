#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2013 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file unittests/libtests/feassemble/data/Mesh1DLinear.odb
##
## @brief Python container holding mesh information for a 1-D
## linear finite-element cell used in testing finite-element C++
## routines.

from pyre.components.Component import Component

import numpy

# ----------------------------------------------------------------------

# Mesh1DLinear class
class Mesh1DLinear(Component):
  """
  Python container holding mesh information for a 1-D linear
  finite-element cell used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="mesh1dlinear"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="mesh")
    
    self.spaceDim = 1
    self.cellDim = 1
    self.numVertices = 2
    self.numCells = 1    
    self.gravityVec = numpy.array( [-1.0e+8], dtype=numpy.float64)
    self.vertices = numpy.array( [[-0.25], [2.0]], dtype=numpy.float64)
    self.cells = numpy.array( [[0, 1]], dtype=numpy.int32)
    self.verticesRef = numpy.array( [[-1.0], [+1.0]], dtype=numpy.float64)

    self.minCellWidth = 2.25
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def mesh():
  """
  Factory for Mesh1DLinear.
  """
  return Mesh1DLinear()


# End of file 
