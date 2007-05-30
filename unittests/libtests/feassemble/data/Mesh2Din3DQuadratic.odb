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

## @file unittests/libtests/feassemble/data/Mesh2Din3DQuadratic.odb
##
## @brief Python container holding mesh information for a 2-D
## quadratic finite-element cell used in testing finite-element C++
## routines.

from pyre.components.Component import Component

import numpy

# ----------------------------------------------------------------------

# Mesh2Din3DQuadratic class
class Mesh2Din3DQuadratic(Component):
  """
  Python container holding mesh information for a 2-D quadratic
  finite-element cell used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="mesh2din3dquadratic"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="mesh")
    
    self.spaceDim = 3
    self.cellDim = 2
    self.numVertices = 6
    self.numCells = 1
    self.vertices = numpy.array( [[ 2.0, -0.5, -0.5],
                                  [ 0.5,  3.0,  0.0],
                                  [-0.5,  0.0,  2.0],
                                  [ 1.3,  1.2, -0.3],
                                  [ 0.1,  1.4,  0.9],
                                  [ 0.8, -0.3,  0.7]],
                                 dtype=numpy.float64)
    self.cells = numpy.array( [[0, 1, 2, 3, 4, 5]], dtype=numpy.int32)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def mesh():
  """
  Factory for Mesh2Din3DQuadratic.
  """
  return Mesh2Din3DQuadratic()


# End of file 
