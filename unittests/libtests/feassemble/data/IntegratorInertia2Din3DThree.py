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

## @file unittests/libtests/feassemble/data/IntegratorInertia2Din3DThree.py

## @brief Python application for generating C++ data files for testing
## C++ IntegratorInertia object with 2-D cell in 3-D space and linear basis
## functions using 2 cells sharing 2 vertices.

from IntegratorInertia import IntegratorInertia

import numpy

# ----------------------------------------------------------------------

# IntegratorInertia2Din3DThree class
class IntegratorInertia2Din3DThree(IntegratorInertia):
  """
  Python application for generating C++ data files for testing C++
  IntegratorInertia object with 2-D cell in 3-D space and linear
  basis functions using 2 cells sharing 2 vertices.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="integratorinertia2din3dthree"):
    """
    Constructor.
    """
    IntegratorInertia.__init__(self, name)

    from Quadrature2Din3DLinearXYZ import Quadrature2Din3DLinearXYZ
    self.quadrature = Quadrature2Din3DLinearXYZ()
    
    self.numVertices = 4
    self.numCells = 2
    self.fiberDim = 3
    
    self.vertices = numpy.array( [[ 0.6, -2.0, -0.4],
                                  [ 3.0,  0.5,  0.0],
                                  [-1.0,  2.0,  4.0],
                                  [-3.4, -0.5,  3.6]],
                                 dtype=numpy.float64)
    self.cells = numpy.array( [[0, 1, 2],
                               [2, 3, 0]], dtype=numpy.int32)
    self.fieldIn = numpy.array( [[ 1.2], [ 0.1], [-0.3],
                                 [ 0.2], [-0.8], [ 1.2],
                                 [ 1.3], [-0.2], [ 1.7],
                                 [ 1.1], [ 1.4], [ 0.9]], dtype=numpy.float64)
    return
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = IntegratorInertia2Din3DThree()
  app.run()


# End of file 
