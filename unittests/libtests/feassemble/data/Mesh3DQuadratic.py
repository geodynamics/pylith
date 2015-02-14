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
    self.gravityVec = numpy.array( [0.0, 0.0, -1.0e8],
                                   dtype=numpy.float64)
    self.vertices = numpy.array( [[-0.5, -2.0, -1.0],
                                  [ 1.0,  1.0,  0.0],
                                  [ 2.0, -2.0, -0.5],
                                  [ 0.2,  0.5,  2.0],
                                  [ 1.5, -0.5, -0.25],
                                  [ 0.25, -0.5, -0.5],
                                  [0.75, -2.0, -0.75],
                                  [-0.15, -0.75,  0.5],
                                  [ 1.1, -0.75,  0.75],
                                  [ 0.6,  0.75,  1.0]],
                                 dtype=numpy.float64)
    self.cells = numpy.array( [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]],
                              dtype=numpy.int32)
    self.verticesRef = numpy.array( [[-1.0, -1.0, -1.0],
                                     [-1.0, +1.0, -1.0],
                                     [+1.0, -1.0, -1.0],
                                     [-1.0, -1.0, +1.0],
                                     [ 0.0,  0.0, -1.0],
                                     [-1.0,  0.0, -1.0],
                                     [ 0.0, -1.0, -1.0],
                                     [-1.0, -1.0,  0.0],
                                     [ 0.0, -1.0,  0.0],
                                     [-1.0,  0.0,  0.0]],
                                    dtype=numpy.float64)

    v0 = self.vertices[0,:]
    v1 = self.vertices[1,:]
    v2 = self.vertices[2,:]
    v3 = self.vertices[3,:]

    e01 = ((v0[0]-v1[0])**2 + (v0[1]-v1[1])**2 + (v0[2]-v1[2])**2)**0.5
    e12 = ((v1[0]-v2[0])**2 + (v1[1]-v2[1])**2 + (v1[2]-v2[2])**2)**0.5
    e20 = ((v2[0]-v0[0])**2 + (v2[1]-v0[1])**2 + (v2[2]-v0[2])**2)**0.5
    e03 = ((v0[0]-v3[0])**2 + (v0[1]-v3[1])**2 + (v0[2]-v3[2])**2)**0.5
    e13 = ((v1[0]-v3[0])**2 + (v1[1]-v3[1])**2 + (v1[2]-v3[2])**2)**0.5
    e23 = ((v2[0]-v3[0])**2 + (v2[1]-v3[1])**2 + (v2[2]-v3[2])**2)**0.5

    vol = 1.0/6.0*numpy.linalg.det(numpy.array([[1.0, v0[0], v0[1], v0[2]],
                                                [1.0, v2[0], v2[1], v2[2]],
                                                [1.0, v1[0], v1[1], v1[2]],
                                                [1.0, v3[0], v3[1], v3[2]]],
                                               dtype=numpy.float64))
    cross012 = numpy.cross(v1-v0, v2-v0)
    area012 = 0.5*(numpy.dot(cross012, cross012))**0.5

    cross013 = numpy.cross(v1-v0, v3-v0)
    area013 = 0.5*(numpy.dot(cross013, cross013))**0.5

    cross123 = numpy.cross(v2-v1, v3-v1)
    area123 = 0.5*(numpy.dot(cross123, cross123))**0.5

    cross203 = numpy.cross(v0-v2, v3-v2)
    area203 = 0.5*(numpy.dot(cross203, cross203))**0.5

    area = area012 + area013 + area123 + area203;
    r = 3.0 * vol / area
    self.minCellWidth = min(e01, e12, e20, e03, e13, e23, 6.38*r)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def mesh():
  """
  Factory for Mesh3DQuadratic.
  """
  return Mesh3DQuadratic()


# End of file 
