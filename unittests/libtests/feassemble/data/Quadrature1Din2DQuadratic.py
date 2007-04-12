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

## @file unittests/libtests/feassemble/data/Quadrature1Din2DQuadratic.py

## @brief Python application for generating C++ data files for testing
## C++ Quadrature1Din2D object w/quadratic basis functions.

from Quadrature1DQuadratic import *

import numpy

# ----------------------------------------------------------------------

# Quadrature1Din2DQuadratic class
class Quadrature1Din2DQuadratic(Quadrature1DQuadratic):
  """
  Python application for generating C++ data files for testing C++
  Quadrature1Din2D object w/quadratic basis functions.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="quadrature1din2dquadratic"):
    """
    Constructor.
    """
    Quadrature1DQuadratic.__init__(self, name)
    
    self.numVertices = 3
    self.spaceDim = 2
    self.numCells = 1
    self.cellDim = 1
    self.numBasis = 3
    self.numQuadPts = 2
    
    self.vertices = numpy.array( [[-0.2, -0.5], [0.3, -0.2], [0.7, 0.3]],
                                 dtype=numpy.float64)
    self.cells = numpy.array( [[0, 1, 2]], dtype=numpy.int32)
    
    return


# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = Quadrature1Din2DQuadratic()
  app.run()


# End of file 
