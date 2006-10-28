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

## @file unittests/libtests/feassemble/data/IntegratorInertia1DLinear.py

## @brief Python application for generating C++ data files for testing
## C++ IntegratorInertia object with 1-D cell and linear basis
## functions.

from IntegratorInertia import IntegratorInertia

import numpy

# ----------------------------------------------------------------------

# IntegratorInertia1DLinear class
class IntegratorInertia1DLinear(IntegratorInertia):
  """
  Python application for generating C++ data files for testing C++
  IntegratorInertia object with 1-D cell and linear basis functions.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="itnegratorinertia1dlinear"):
    """
    Constructor.
    """
    IntegratorInertia.__init__(self, name)

    from Quadrature1DLinear import Quadrature1DLinear
    self.quadrature = Quadrature1DLinear()
    
    self.numVertices = 2
    self.numCells = 1
    self.fiberDim = 1
    
    self.vertices = numpy.array( [[-0.25], [2.0]], dtype=numpy.float64)
    self.cells = numpy.array( [[0, 1]], dtype=numpy.int32)
    self.fieldIn = numpy.array( [[1.2], [1.5]], dtype=numpy.float64)
    return
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = IntegratorInertia1DLinear()
  app.run()


# End of file 
