#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file unittests/libtests/feassemble/data/Mesh2Din3DLinearYZ.odb
##
## @brief Python container holding mesh information for a 2-D
## linear finite-element cell used in testing finite-element C++
## routines.

from pyre.components.Component import Component

import numpy

# ----------------------------------------------------------------------

# Mesh2Din3DLinearYZ class
class Mesh2Din3DLinearYZ(Component):
  """
  Python container holding mesh information for a 2-D linear
  finite-element cell used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="mesh2din3dlinearyz"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="mesh")
    
    self.spaceDim = 3
    self.cellDim = 2
    self.numVertices = 3
    self.numCells = 1
    self.gravityVec = numpy.array( [0.0, 0.0, -1.0e8],
                                   dtype=numpy.float64)
    self.vertices = numpy.array( [[ 0.0,  0.0,  0.0],
                                  [ 0.0,  1.0,  0.0],
                                  [ 0.0,  0.0,  1.0]],
                                 dtype=numpy.float64)
    self.cells = numpy.array( [[0, 1, 2]], dtype=numpy.int32)
    self.verticesRef = numpy.array( [[-1.0, -1.0],
                                     [+1.0, -1.0],
                                     [-1.0, +1.0]],
                                    dtype=numpy.float64)

    a = (0.0**2 + 1.0**2 + 0.0**2)**0.5
    b = (0.0**2 + 1.0**2 + 0.0**2)**0.5
    c = (0.0**2 + 0.0**2 + 1.0**2)**0.5
    k = 0.5 * (a + b + c)
    r = (k*(k-a)*(k-b)*(k-c))**0.5 / k
    self.minCellWidth = r
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def mesh():
  """
  Factory for Mesh2Din3DLinearYZ.
  """
  return Mesh2Din3DLinearYZ()


# End of file 
