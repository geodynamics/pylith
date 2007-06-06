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

## @file unittests/libtests/feassemble/data/Mesh1Din2DQuadratic.odb
##
## @brief Python container holding mesh information for a 1-D
## quadratic finite-element cell used in testing finite-element C++
## routines.

from pyre.components.Component import Component

import numpy

# ----------------------------------------------------------------------

# Mesh1Din2DQuadratic class
class Mesh1Din2DQuadratic(Component):
  """
  Python container holding mesh information for a 1-D quadratic
  finite-element cell used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="mesh1din2dquadratic"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="mesh")
    
    self.spaceDim = 2
    self.cellDim = 1
    self.numVertices = 3
    self.numCells = 1
    self.vertices = numpy.array( [[-0.2, -0.5], [0.3, -0.2], [0.7, 0.3]],
                                 dtype=numpy.float64)
    self.cells = numpy.array( [[0, 1, 2]], dtype=numpy.int32)
    self.verticesRef = numpy.array([[-1.0], [1.0], [0.0]], dtype=numpy.float64)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def mesh():
  """
  Factory for Mesh1Din2DQuadratic.
  """
  return Mesh1Din2DQuadratic()


# End of file 
