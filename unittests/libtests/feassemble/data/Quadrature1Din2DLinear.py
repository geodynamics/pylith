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

## @file unittests/libtests/feassemble/data/Quadrature1Din2DLinear.py

## @brief Python application for generating C++ data files for testing
## C++ Quadrature1Din2D object w/linear basis functions.

from Quadrature1DLinear import *

import numpy

# ----------------------------------------------------------------------

# Quadrature1Din2DLinear class
class Quadrature1Din2DLinear(Quadrature1DLinear):
  """
  Python application for generating C++ data files for testing C++
  Quadrature1Din2D object w/linear basis functions.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="quadrature1din2dlinear"):
    """
    Constructor.
    """
    Quadrature1DLinear.__init__(self, name)
    
    self.numVertices = 2
    self.spaceDim = 2
    self.numCells = 1
    self.cellDim = 1
    self.numBasis = 2
    self.numQuadPts = 1
    
    self.vertices = numpy.array( [[-0.2, -0.5], [0.7, 0.3]],
                                 dtype=numpy.float64)
    self.cells = numpy.array( [[0, 1]], dtype=numpy.int32)
    return
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = Quadrature1Din2DLinear()
  app.run()


# End of file 
