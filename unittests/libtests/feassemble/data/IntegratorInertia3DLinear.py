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

## @file unittests/libtests/feassemble/data/IntegratorInertia3DLinear.py

## @brief Python application for generating C++ data files for testing
## C++ IntegratorInertia object with 3-D cell and linear basis
## functions.

from IntegratorInertia import IntegratorInertia

import numpy

# ----------------------------------------------------------------------

# IntegratorInertia3DLinear class
class IntegratorInertia3DLinear(IntegratorInertia):
  """
  Python application for generating C++ data files for testing C++
  IntegratorInertia object with 3-D cell and linear basis functions.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="itnegratorinertia3dlinear"):
    """
    Constructor.
    """
    IntegratorInertia.__init__(self, name)

    from Quadrature3DLinear import Quadrature3DLinear
    self.quadrature = Quadrature3DLinear()
    
    self.numVertices = 4
    self.numCells = 1
    self.fiberDim = 3
    
    self.vertices = numpy.array( [[-0.5, -1.0, -0.5],
                                  [ 2.0, -0.5, -0.4],
                                  [ 1.0, -0.1, -0.3],
                                  [-0.2,  0.5,  2.0]], dtype=numpy.float64)
    self.cells = numpy.array( [[0, 1, 2, 3]], dtype=numpy.Int32)
    self.fieldIn = numpy.array( [[ 1.2], [ 0.1], [-0.3],
                                 [ 0.2], [-0.8], [ 1.2],
                                 [ 1.3], [-0.2], [ 1.7],
                                 [ 1.1], [ 1.4], [ 0.9]], dtype=numpy.float64)
    return
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = IntegratorInertia3DLinear()
  app.run()


# End of file 
