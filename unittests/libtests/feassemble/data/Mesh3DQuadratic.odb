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

## @file unittests/libtests/feassemble/data/Mesh3DQuadratic.odb
##
## @brief Python container holding mesh information for a 2-D
## quadratic finite-element cell used in testing finite-element C++
## routines.

from pyre.components.Component import Component

import numpy

# ----------------------------------------------------------------------

# Mesh3DQuadratic class
class Mesh3DQuadratic(Component):
  """
  Python container holding mesh information for a 1-D quadratic
  finite-element cell used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="mesh3dquadratic"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="mesh")
    
    self.spaceDim = 3
    self.cellDim = 3
    self.numVertices = 10
    self.numCells = 1
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
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def mesh():
  """
  Factory for Mesh3DQuadratic.
  """
  return Mesh3DQuadratic()


# End of file 
