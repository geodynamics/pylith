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

## @file unittests/libtests/feassemble/data/IntegratorInertia3DQuadratic.py

## @brief Python application for generating C++ data files for testing
## C++ IntegratorInertia object with 1-D cell and quadratic basis
## functions.

from IntegratorInertia import IntegratorInertia

import numpy

# ----------------------------------------------------------------------

# IntegratorInertia3DQuadratic class
class IntegratorInertia3DQuadratic(IntegratorInertia):
  """
  Python application for generating C++ data files for testing C++
  IntegratorInertia object with 1-D cell and quadratic basis functions.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="integratorinertia3dquadratic"):
    """
    Constructor.
    """
    IntegratorInertia.__init__(self, name)

    from Quadrature3DQuadratic import Quadrature3DQuadratic
    self.quadrature = Quadrature3DQuadratic()
    
    self.numVertices = 10
    self.numCells = 1
    self.fiberDim = 3
    
    self.vertices = numpy.array( [[-0.5, -2.0, -1.0],
                                  [ 2.0, -2.0, -0.5],
                                  [ 1.0,  1.0,  0.0],
                                  [ 0.2,  0.5,  2.0],
                                  [ 0.7, -2.1, -0.8],
                                  [ 0.3, -0.5, -0.5],
                                  [-0.2, -0.8,  0.5],
                                  [ 1.5, -0.6, -0.2],
                                  [ 0.6,  0.8,  0.9],
                                  [ 1.1, -0.8,  0.7]],
                                 dtype=numpy.float64)
    self.cells = numpy.array( [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]],
                              dtype=numpy.int32)
    self.fieldIn = numpy.array( [[1.2], [1.5], [0.3],
                                 [0.8], [0.9], [1.4],
                                 [0.7], [0.2], [0.8],
                                 [1.4], [1.5], [1.4],
                                 [0.5], [0.7], [1.6],
                                 [1.9], [1.2], [1.3],
                                 [0.6], [0.3], [1.5],
                                 [1.3], [0.6], [0.3],
                                 [1.4], [0.1], [1.8],
                                 [0.2], [1.0], [0.8]], dtype=numpy.float64)
    return
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = IntegratorInertia3DQuadratic()
  app.run()


# End of file 
