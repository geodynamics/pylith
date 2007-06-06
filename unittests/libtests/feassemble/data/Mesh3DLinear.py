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

## @file unittests/libtests/feassemble/data/Mesh3DLinear.odb
##
## @brief Python container holding mesh information for a 3-D
## linear finite-element cell used in testing finite-element C++
## routines.

from pyre.components.Component import Component

import numpy

# ----------------------------------------------------------------------

# Mesh3DLinear class
class Mesh3DLinear(Component):
  """
  Python container holding mesh information for a 3-D linear
  finite-element cell used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="mesh3dlinear"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="mesh")
    
    self.spaceDim = 3
    self.cellDim = 3
    self.numVertices = 4
    self.numCells = 1
    self.vertices = numpy.array( [[-0.5, -1.0, -0.5],
                                  [ 2.0, -0.5, -0.4],
                                  [ 1.0, -0.1, -0.3],
                                  [-0.2,  0.5,  2.0]], dtype=numpy.float64)
    self.cells = numpy.array( [[0, 1, 2, 3]], dtype=numpy.int32)
    self.verticesRef = numpy.array( [[-1.0, -1.0, -1.0],
                                     [+1.0, -1.0, -1.0],
                                     [-1.0, +1.0, -1.0],
                                     [-1.0, -1.0, +1.0]],
                                    dtype=numpy.float64)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def mesh():
  """
  Factory for Mesh3DLinear.
  """
  return Mesh3DLinear()


# End of file 
