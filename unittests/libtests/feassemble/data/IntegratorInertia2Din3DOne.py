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

## @file unittests/libtests/feassemble/data/IntegratorInertia2Din3DOne.py

## @brief Python application for generating C++ data files for testing
## C++ IntegratorInertia object with 2-D cell in 3-D space and linear basis
## functions.

from IntegratorInertia import IntegratorInertia

import numpy

# ----------------------------------------------------------------------

# IntegratorInertia2Din3DOne class
class IntegratorInertia2Din3DOne(IntegratorInertia):
  """
  Python application for generating C++ data files for testing C++
  IntegratorInertia object with 2-D cell in 3-D space and linear
  basis functions.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="integratorinertia2din3done"):
    """
    Constructor.
    """
    IntegratorInertia.__init__(self, name)

    from Quadrature2Din3DLinearXYZ import Quadrature2Din3DLinearXYZ
    self.quadrature = Quadrature2Din3DLinearXYZ()
    
    self.numVertices = 3
    self.numCells = 1
    self.fiberDim = 3
    
    self.vertices = numpy.array( [[ 0.5, -2.0, -0.5],
                                  [ 3.0,  0.5,  0.0],
                                  [-1.0,  2.0,  4.0]],
                                 dtype=numpy.float64)
    self.cells = numpy.array( [[0, 1, 2]], dtype=numpy.int32)
    self.fieldIn = numpy.array( [[1.2], [0.1], [-0.3],
                                 [0.5], [-0.3], [1.2],
                                 [1.1], [0.5], [0.8]], dtype=numpy.float64)
    return
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = IntegratorInertia2Din3DOne()
  app.run()


# End of file 
