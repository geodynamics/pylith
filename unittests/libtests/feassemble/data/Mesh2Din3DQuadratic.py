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
    self.gravityVec = numpy.array( [0.0, 0.0, -1.0e8],
                                   dtype=numpy.float64)
    self.vertices = numpy.array( [[ 2.0, -0.5, -0.5],
                                  [ 0.5,  3.0,  0.0],
                                  [-0.5,  0.0,  2.0],
                                  [ 0.0,  1.5,  1.0],
                                  [ 0.75, -0.25, 0.75],
                                  [ 1.25, 1.25, -0.25]],
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
  Factory for Mesh2Din3DQuadratic.
  """
  return Mesh2Din3DQuadratic()


# End of file 
