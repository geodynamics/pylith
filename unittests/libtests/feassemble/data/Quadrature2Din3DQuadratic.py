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

## @file unittests/libtests/feassemble/data/Quadrature2Din3DQuadratic.py

## @brief Python application for generating C++ data files for testing
## C++ Quadrature2Din3D object w/quadratic basis functions.

from Quadrature2DQuadratic import *

import numpy

# ----------------------------------------------------------------------

# Quadrature2Din3DQuadratic class
class Quadrature2Din3DQuadratic(Quadrature2DQuadratic):
  """
  Python application for generating C++ data files for testing C++
  Quadrature2Din3D object w/quadratic basis functions.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="quadrature2din3dquadratic"):
    """
    Constructor.
    """
    Quadrature2DQuadratic.__init__(self, name)
    
    self.numVertices = 6
    self.spaceDim = 3
    self.numCells = 1
    self.cellDim = 2
    self.numCorners = 6
    self.numQuadPts = 3
    
    self.vertices = numpy.array( [[ 2.0, -0.5, -0.5],
                                  [ 0.5,  3.0,  0.0],
                                  [-0.5,  0.0,  2.0],
                                  [ 1.3,  1.2, -0.3],
                                  [ 0.1,  1.4,  0.9],
                                  [ 0.8, -0.3,  0.7]],
                                 dtype=numpy.Float64)
    self.cells = numpy.array( [[0, 1, 2, 3, 4, 5]], dtype=numpy.Int32)
    
    return


# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = Quadrature2Din3DQuadratic()
  app.run()


# End of file 
