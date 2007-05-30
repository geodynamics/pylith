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

## @file unittests/libtests/feassemble/data/Mesh1Din2DLinear.odb
##
## @brief Python container holding mesh information for a 1-D
## linear finite-element cell used in testing finite-element C++
## routines.

from pyre.components.Component import Component

import numpy

# ----------------------------------------------------------------------

# Mesh1Din2DLinear class
class Mesh1Din2DLinear(Component):
  """
  Python container holding mesh information for a 1-D linear
  finite-element cell used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="mesh1din2dlinear"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="mesh")
    
    self.spaceDim = 2
    self.cellDim = 1
    self.numVertices = 2
    self.numCells = 1
    self.vertices = numpy.array( [[-0.2, -0.5], [0.7, 0.3]],
                                 dtype=numpy.float64)
    self.cells = numpy.array( [[0, 1]], dtype=numpy.int32)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def mesh():
  """
  Factory for Mesh1Din2DLinear.
  """
  return Mesh1Din2DLinear()


# End of file 
