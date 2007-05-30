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

## @file unittests/libtests/feassemble/data/Mesh2Din3DLinearXY.odb
##
## @brief Python container holding mesh information for a 2-D
## linear finite-element cell used in testing finite-element C++
## routines.

from pyre.components.Component import Component

import numpy

# ----------------------------------------------------------------------

# Mesh2Din3DLinearXY class
class Mesh2Din3DLinearXY(Component):
  """
  Python container holding mesh information for a 2-D linear
  finite-element cell used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="mesh2din3dlinearxy"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="mesh")
    
    self.spaceDim = 3
    self.cellDim = 2
    self.numVertices = 3
    self.numCells = 1
    self.vertices = numpy.array( [[ 0.0,  0.0,  0.0],
                                  [-1.0,  0.0,  0.0],
                                  [ 0.0,  0.0,  1.0]],
                                 dtype=numpy.float64)
    self.cells = numpy.array( [[0, 1, 2]], dtype=numpy.int32)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def mesh():
  """
  Factory for Mesh2Din3DLinearXY.
  """
  return Mesh2Din3DLinearXY()


# End of file 
