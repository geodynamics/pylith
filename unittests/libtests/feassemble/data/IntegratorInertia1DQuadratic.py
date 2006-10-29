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

## @file unittests/libtests/feassemble/data/IntegratorInertia1DQuadratic.py

## @brief Python application for generating C++ data files for testing
## C++ IntegratorInertia object with 1-D cell and quadratic basis
## functions.

from IntegratorInertia import IntegratorInertia

import numpy

# ----------------------------------------------------------------------

# IntegratorInertia1DQuadratic class
class IntegratorInertia1DQuadratic(IntegratorInertia):
  """
  Python application for generating C++ data files for testing C++
  IntegratorInertia object with 1-D cell and quadratic basis functions.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="integratorinertia1dquadratic"):
    """
    Constructor.
    """
    IntegratorInertia.__init__(self, name)

    from Quadrature1DQuadratic import Quadrature1DQuadratic
    self.quadrature = Quadrature1DQuadratic()
    
    self.numVertices = 3
    self.numCells = 1
    self.fiberDim = 1
    
    self.vertices = numpy.array( [[-0.25], [0.875], [2.0]],
                                 dtype=numpy.float64)
    self.cells = numpy.array( [[0, 1, 2]], dtype=numpy.Int32)
    self.fieldIn = numpy.array( [[1.2], [1.5], [-0.8]], dtype=numpy.float64)
    return
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = IntegratorInertia1DQuadratic()
  app.run()


# End of file 
