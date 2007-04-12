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

## @file unittests/libtests/feassemble/data/Quadrature1Din3DQuadratic.py

## @brief Python application for generating C++ data files for testing
## C++ Quadrature1Din3D object w/quadratic basis functions.

from Quadrature1DQuadratic import *

import numpy

# ----------------------------------------------------------------------

# Quadrature1Din3DQuadratic class
class Quadrature1Din3DQuadratic(Quadrature1DQuadratic):
  """
  Python application for generating C++ data files for testing C++
  Quadrature1Din3D object w/quadratic basis functions.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="quadrature1din3dquadratic"):
    """
    Constructor.
    """
    Quadrature1DQuadratic.__init__(self, name)
    
    self.numVertices = 3
    self.spaceDim = 3
    self.numCells = 1
    self.cellDim = 1
    self.numBasis = 3
    self.numQuadPts = 2
    
    self.vertices = numpy.array( [[1.0, -1.5, -2.0],
                                  [0.3, 0.3, 0.8],
                                  [-0.5, 2.0, 3.0]],
                                 dtype=numpy.float64)
    self.cells = numpy.array( [[0, 1, 2]], dtype=numpy.int32)
    
    return


# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = Quadrature1Din3DQuadratic()
  app.run()


# End of file 
