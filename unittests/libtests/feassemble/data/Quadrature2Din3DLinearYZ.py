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

## @file unittests/libtests/feassemble/data/Quadrature2Din3DLinearYZ.py

## @brief Python application for generating C++ data files for testing
## C++ Quadrature2Din3D object w/linear basis functions.

from Quadrature2DLinear import *

import numpy

# ----------------------------------------------------------------------

# Quadrature2Din3DLinearYZ class
class Quadrature2Din3DLinearYZ(Quadrature2DLinear):
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
    
    self.vertices = numpy.array( [[ 0.0,  0.0,  0.0],
                                  [ 0.0,  1.0,  0.0],
                                  [ 0.0,  0.0,  1.0]],
                                 dtype=numpy.Float64)
    self.cells = numpy.array( [[0, 1, 2]], dtype=numpy.Int32)
    return
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = Quadrature2Din3DLinearYZ()
  app.run()


# End of file 
