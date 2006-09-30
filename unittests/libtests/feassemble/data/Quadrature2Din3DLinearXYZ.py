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

## @file unittests/libtests/feassemble/data/Quadrature2Din3DLinearXYZ.py

## @brief Python application for generating C++ data files for testing
## C++ Quadrature2Din3D object w/linear basis functions.

from Quadrature2DLinear import *

import numpy

# ----------------------------------------------------------------------

# Quadrature2Din3DLinearXYZ class
class Quadrature2Din3DLinearXYZ(Quadrature2DLinear):
  """
  Python application for generating C++ data files for testing C++
  Quadrature2Din3D object w/linear basis functions.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="quadrature2din3dlinear"):
    """
    Constructor.
    """
    Quadrature2DLinear.__init__(self, name)
    
    self.numVertices = 3
    self.spaceDim = 3
    self.numCells = 1
    self.cellDim = 2
    self.numCorners = 3
    self.numQuadPts = 1
    
    self.vertices = numpy.array( [[ 0.5, -2.0, -0.5],
                                  [ 3.0,  0.5,  0.0],
                                  [-1.0,  2.0,  4.0]],
                                 dtype=numpy.Float64)
    self.cells = numpy.array( [[0, 1, 2]], dtype=numpy.Int32)
    return
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = Quadrature2Din3DLinearXYZ()
  app.run()


# End of file 
