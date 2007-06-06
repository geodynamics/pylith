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

## @file unittests/libtests/feassemble/data/Mesh2DQuadratic.odb
##
## @brief Python container holding mesh information for a 2-D
## quadratic finite-element cell used in testing finite-element C++
## routines.

from pyre.components.Component import Component

import numpy

# ----------------------------------------------------------------------

# Mesh2DQuadratic class
class Mesh2DQuadratic(Component):
  """
  Python container holding mesh information for a 2-D quadratic
  finite-element cell used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="mesh2dquadratic"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="mesh")
    
    self.spaceDim = 2
    self.cellDim = 2
    self.numVertices = 6
    self.numCells = 1
    self.vertices = numpy.array( [[-1.0, -1.0],
                                  [+1.0, +0.2],
                                  [-1.5, +0.5],
                                  [+0.0, -0.6],
                                  [+0.25, +0.35],
                                  [-1.25, -0.25]],
                                 dtype=numpy.float64)
    self.cells = numpy.array( [[0, 1, 2, 3, 4, 5]], dtype=numpy.int32)
    self.verticesRef = numpy.array( [[-1.0, -1.0],
                                     [+1.0, -1.0],
                                     [-1.0, +1.0],
                                     [ 0.0,  0.0],
                                     [-1.0,  0.0],
                                     [ 0.0, -1.0]],
                                    dtype=numpy.float64)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def mesh():
  """
  Factory for Mesh2DQuadratic.
  """
  return Mesh2DQuadratic()


# End of file 
