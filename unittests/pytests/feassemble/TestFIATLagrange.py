#!/usr/bin/env python
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file unittests/pytests/feassemble/TestFIATLagrange.py

## @brief Unit testing of FIATLagrange object.

import unittest
import numpy
from pylith.feassemble.FIATLagrange import FIATLagrange
from pylith.utils.testarray import test_double

# ----------------------------------------------------------------------
class Line2(object):

  def __init__(self):
    """
    Setup line2 cell.
    """
    vertices = numpy.array([[-1.0], [1.0]])
    quadPts = numpy.array( [[0.0]],
                                dtype=numpy.float64 )
    quadWts = numpy.array( [2.0], dtype=numpy.float64 )

    # Compute basis fns and derivatives at quadrature points
    basis = numpy.zeros( (1, 2), dtype=numpy.float64)
    basisDeriv = numpy.zeros( (1, 2, 1), dtype=numpy.float64)
    iQuad = 0
    for q in quadPts:
      basis[iQuad] = numpy.array([self.N0(q), self.N1(q)],
                                 dtype=numpy.float64).reshape( (2,) )
      deriv = numpy.array([[self.N0p(q)], [self.N1p(q)]],
                          dtype=numpy.float64)      
      basisDeriv[iQuad] = deriv.reshape((2, 1))
      iQuad += 1

    self.cellDim = 1
    self.numCorners = len(vertices)
    self.numQuadPts = len(quadPts)
    self.vertices = vertices
    self.quadPts = quadPts
    self.quadWts = quadWts
    self.basis = basis
    self.basisDeriv = basisDeriv
    return


  def N0(self, p):
    return 0.5*(1.0-p)

  def N0p(self, p):
    return -0.5

  def N1(self, p):
    return 0.5*(1.0+p)

  def N1p(self, p):
    return 0.5


# ----------------------------------------------------------------------
class Line3(object):

  def __init__(self):
    """
    Setup line3 cell.
    """
    vertices = numpy.array([[-1.0], [1.0], [0.0]])
    quadPts = numpy.array([ [-1.0/3**0.5],
                            [+1.0/3**0.5] ])
    quadWts = numpy.array( [1.0, 1.0])

    # Compute basis fns and derivatives at quadrature points
    basis = numpy.zeros( (2, 3), dtype=numpy.float64)
    basisDeriv = numpy.zeros( (2, 3, 1), dtype=numpy.float64)
    iQuad = 0
    for q in quadPts:
      basis[iQuad] = numpy.array([self.N0(q), self.N1(q), self.N2(q)],
                                 dtype=numpy.float64).reshape( (3,) )
      deriv = numpy.array([[self.N0p(q)], [self.N1p(q)], [self.N2p(q)]])      
      basisDeriv[iQuad] = deriv.reshape((3, 1))
      iQuad += 1

    self.cellDim = 1
    self.numCorners = len(vertices)
    self.numQuadPts = len(quadPts)
    self.vertices = vertices
    self.quadPts = quadPts
    self.quadWts = quadWts
    self.basis = basis
    self.basisDeriv = basisDeriv
    return


  def N0(self, p):
    return -0.5*p*(1.0-p)

  def N0p(self, p):
    return 1.0*p - 0.5

  def N1(self, p):
    return 0.5*p*(1.0+p)

  def N1p(self, p):
    return 1.0*p + 0.5

  def N2(self, p):
    return (1.0-p*p)

  def N2p(self, p):
    return -2.0*p


# ----------------------------------------------------------------------
class Quad4(object):

  def __init__(self):
    """
    Setup quad4 cell.
    """
    vertices = numpy.array([[-1.0, -1.0],
                            [+1.0, -1.0],
                            [+1.0, +1.0],
                            [-1.0, +1.0]])
    quadPts = numpy.array([ [-1.0/3**0.5, -1.0/3**0.5],
                            [+1.0/3**0.5, -1.0/3**0.5],
                            [-1.0/3**0.5, +1.0/3**0.5],
                            [+1.0/3**0.5, +1.0/3**0.5] ])
    quadWts = numpy.array( [1.0, 1.0, 1.0, 1.0])

    # Compute basis fns and derivatives at quadrature points
    basis = numpy.zeros( (4, 4), dtype=numpy.float64)
    basisDeriv = numpy.zeros( (4, 4, 2), dtype=numpy.float64)
    iQuad = 0
    for q in quadPts:
      basis[iQuad] = numpy.array([self.N0(q), self.N1(q),
                                  self.N2(q), self.N3(q)],
                                 dtype=numpy.float64).reshape( (4,) )
      deriv = numpy.array([[self.N0p(q), self.N0q(q)],
                           [self.N1p(q), self.N1q(q)],
                           [self.N2p(q), self.N2q(q)],
                           [self.N3p(q), self.N3q(q)]])
      basisDeriv[iQuad] = deriv.reshape((4, 2))
      iQuad += 1

    self.cellDim = 2
    self.numCorners = len(vertices)
    self.numQuadPts = len(quadPts)
    self.vertices = vertices
    self.quadPts = quadPts
    self.quadWts = quadWts
    self.basis = basis
    self.basisDeriv = basisDeriv
    return


  def N0(self, p):
    return (1-p[0])*(1-p[1])/4.0

  def N0p(self, p):
    return -(1-p[1])/4.0

  def N0q(self, p):
    return -(1-p[0])/4.0

  def N1(self, p):
    return (1+p[0])*(1-p[1])/4.0

  def N1p(self, p):
    return (1-p[1])/4.0

  def N1q(self, p):
    return -(1+p[0])/4.0

  def N2(self, p):
    return (1+p[0])*(1+p[1])/4.0

  def N2p(self, p):
    return (1+p[1])/4.0

  def N2q(self, p):
    return (1+p[0])/4.0

  def N3(self, p):
    return (1-p[0])*(1+p[1])/4.0

  def N3p(self, p):
    return -(1+p[1])/4.0

  def N3q(self, p):
    return (1-p[0])/4.0


# ----------------------------------------------------------------------
class Quad9(object):

  def __init__(self):
    """
    Setup quad9 cell.
    """
    vertices = numpy.array([[-1.0, -1.0], # corners
                            [+1.0, -1.0],
                            [+1.0, +1.0],
                            [-1.0, +1.0],
                            [ 0.0, -1.0], # edges
                            [+1.0,  0.0],
                            [ 0.0, +1.0],
                            [-1.0,  0.0],
                            [ 0.0,  0.0], # face
                            ])
    quadPts = numpy.array([ [-(3.0/5)**0.5, -(3.0/5)**0.5],
                            [          0.0, -(3.0/5)**0.5],
                            [+(3.0/5)**0.5, -(3.0/5)**0.5],
                            [-(3.0/5)**0.5, 0.0],
                            [          0.0, 0.0],
                            [+(3.0/5)**0.5, 0.0],
                            [-(3.0/5)**0.5, +(3.0/5)**0.5],
                            [          0.0, +(3.0/5)**0.5],
                            [+(3.0/5)**0.5, +(3.0/5)**0.5] ])
    quadWts = numpy.array( [25.0/81, 40.0/81, 25.0/81,
                            40.0/81, 64.0/81, 40.0/81,
                            25.0/81, 40.0/81, 25.0/81])

    # Compute basis fns and derivatives at quadrature points
    basis = numpy.zeros( (9, 9), dtype=numpy.float64)
    basisDeriv = numpy.zeros( (9, 9, 2), dtype=numpy.float64)
    iQuad = 0
    for q in quadPts:
      basis[iQuad] = numpy.array([self.N0(q), self.N1(q), self.N2(q),
                                  self.N3(q), self.N4(q), self.N5(q),
                                  self.N6(q), self.N7(q), self.N8(q)],
                                 dtype=numpy.float64).reshape( (9,) )
      print basis[iQuad]
      deriv = numpy.array([[self.N0p(q), self.N0q(q)],
                           [self.N1p(q), self.N1q(q)],
                           [self.N2p(q), self.N2q(q)],
                           [self.N3p(q), self.N3q(q)],
                           [self.N4p(q), self.N4q(q)],
                           [self.N5p(q), self.N5q(q)],
                           [self.N6p(q), self.N6q(q)],
                           [self.N7p(q), self.N7q(q)],
                           [self.N8p(q), self.N8q(q)]])
      basisDeriv[iQuad] = deriv.reshape((9, 2))
      iQuad += 1

    self.cellDim = 2
    self.numCorners = len(vertices)
    self.numQuadPts = len(quadPts)
    self.vertices = vertices
    self.quadPts = quadPts
    self.quadWts = quadWts
    self.basis = basis
    self.basisDeriv = basisDeriv
    return


  def N0(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(-0.5*p[1])*(1.0-p[1])

  def N0p(self, p):
    return (p[0]-0.5)*(-0.5*p[1])*(1.0-p[1])

  def N0q(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(p[1]-0.5)

  def N1(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(-0.5*p[1])*(1.0-p[1])

  def N1p(self, p):
    return (p[0]+0.5)*(-0.5*p[1])*(1.0-p[1])

  def N1q(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(p[1]-0.5)

  def N2(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(+0.5*p[1])*(1.0+p[1])

  def N2p(self, p):
    return (p[0]+0.5)*(+0.5*p[1])*(1.0+p[1])

  def N2q(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(p[1]+0.5)

  def N3(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(+0.5*p[1])*(1.0+p[1])

  def N3p(self, p):
    return (p[0]-0.5)*(+0.5*p[1])*(1.0+p[1])

  def N3q(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(p[1]+0.5)

  def N4(self, p):
    return (1.0-p[0]**2)*(-0.5*p[1])*(1.0-p[1])

  def N4p(self, p):
    return (-2.0*p[0])*(-0.5*p[1])*(1.0-p[1])

  def N4q(self, p):
    return (1.0-p[0]**2)*(p[1]-0.5)

  def N5(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(1.0-p[1]**2)

  def N5p(self, p):
    return (p[0]+0.5)*(1.0-p[1]**2)

  def N5q(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(-2.0*p[1])

  def N6(self, p):
    return (1.0-p[0]**2)*(+0.5*p[1])*(1.0+p[1])

  def N6p(self, p):
    return (-2.0*p[0])*(+0.5*p[1])*(1.0+p[1])

  def N6q(self, p):
    return (1.0-p[0]**2)*(p[1]+0.5)

  def N7(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(1.0-p[1]**2)

  def N7p(self, p):
    return (p[0]-0.5)*(1.0-p[1]**2)

  def N7q(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(-2.0*p[1])

  def N8(self, p):
    return (1.0-p[0]**2)*(1.0-p[1]**2)

  def N8p(self, p):
    return (-2.0*p[0])*(1.0-p[1]**2)

  def N8q(self, p):
    return (1.0-p[0]**2)*(-2.0*p[1])


# ----------------------------------------------------------------------
class Hex8(object):

  def __init__(self):
    """
    Setup hex8 cell.
    """
    vertices = numpy.array([[-1.0, -1.0, -1.0],
                            [+1.0, -1.0, -1.0],
                            [+1.0, +1.0, -1.0],
                            [-1.0, +1.0, -1.0],
                            [-1.0, -1.0, +1.0],
                            [+1.0, -1.0, +1.0],
                            [+1.0, +1.0, +1.0],
                            [-1.0, +1.0, +1.0]])
    quadPts = numpy.array([ [-1.0/3**0.5, -1.0/3**0.5, -1.0/3**0.5],
                            [+1.0/3**0.5, -1.0/3**0.5, -1.0/3**0.5],
                            [+1.0/3**0.5, +1.0/3**0.5, -1.0/3**0.5],
                            [-1.0/3**0.5, +1.0/3**0.5, -1.0/3**0.5],
                            [-1.0/3**0.5, -1.0/3**0.5, +1.0/3**0.5],
                            [+1.0/3**0.5, -1.0/3**0.5, +1.0/3**0.5],
                            [+1.0/3**0.5, +1.0/3**0.5, +1.0/3**0.5],
                            [-1.0/3**0.5, +1.0/3**0.5, +1.0/3**0.5]])
    quadWts = numpy.array( [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    # Compute basis fns and derivatives at quadrature points
    basis = numpy.zeros( (8, 8), dtype=numpy.float64)
    basisDeriv = numpy.zeros( (8, 8, 3), dtype=numpy.float64)
    iQuad = 0
    for q in quadPts:
      basis[iQuad] = numpy.array([self.N0(q), self.N1(q),
                                  self.N2(q), self.N3(q),
                                  self.N4(q), self.N5(q),
                                  self.N6(q), self.N7(q)],
                                 dtype=numpy.float64).reshape( (8,) )
      deriv = numpy.array([[self.N0p(q), self.N0q(q), self.N0r(q)],
                           [self.N1p(q), self.N1q(q), self.N1r(q)],
                           [self.N2p(q), self.N2q(q), self.N2r(q)],
                           [self.N3p(q), self.N3q(q), self.N3r(q)],
                           [self.N4p(q), self.N4q(q), self.N4r(q)],
                           [self.N5p(q), self.N5q(q), self.N5r(q)],
                           [self.N6p(q), self.N6q(q), self.N6r(q)],
                           [self.N7p(q), self.N7q(q), self.N7r(q)]])      
      basisDeriv[iQuad] = deriv.reshape((8, 3))
      iQuad += 1

    self.cellDim = 3
    self.numCorners = len(vertices)
    self.numQuadPts = len(quadPts)
    self.vertices = vertices
    self.quadPts = quadPts
    self.quadWts = quadWts
    self.basis = basis
    self.basisDeriv = basisDeriv
    return


  def N0(self, p):
    return (1-p[0])*(1-p[1])*(1-p[2])/8.0
  
  def N0p(self, p):
    return -(1-p[1])*(1-p[2])/8.0
  
  def N0q(self, p):
    return -(1-p[0])*(1-p[2])/8.0
  
  def N0r(self, p):
    return -(1-p[0])*(1-p[1])/8.0
  
  def N1(self, p):
    return (1+p[0])*(1-p[1])*(1-p[2])/8.0
  
  def N1p(self, p):
    return (1-p[1])*(1-p[2])/8.0
  
  def N1q(self, p):
    return -(1+p[0])*(1-p[2])/8.0
  
  def N1r(self, p):
    return -(1+p[0])*(1-p[1])/8.0
  
  def N2(self, p):
    return (1+p[0])*(1+p[1])*(1-p[2])/8.0
  
  def N2p(self, p):
    return (1+p[1])*(1-p[2])/8.0
  
  def N2q(self, p):
    return (1+p[0])*(1-p[2])/8.0

  def N2r(self, p):
    return -(1+p[0])*(1+p[1])/8.0

  def N3(self, p):
    return (1-p[0])*(1+p[1])*(1-p[2])/8.0
  
  def N3p(self, p):
    return -(1+p[1])*(1-p[2])/8.0
  
  def N3q(self, p):
    return (1-p[0])*(1-p[2])/8.0
  
  def N3r(self, p):
    return -(1-p[0])*(1+p[1])/8.0
  
  def N4(self, p):
    return (1-p[0])*(1-p[1])*(1+p[2])/8.0

  def N4p(self, p):
    return -(1-p[1])*(1+p[2])/8.0
  
  def N4q(self, p):
    return -(1-p[0])*(1+p[2])/8.0
  
  def N4r(self, p):
    return (1-p[0])*(1-p[1])/8.0
  
  def N5(self, p):
    return (1+p[0])*(1-p[1])*(1+p[2])/8.0
  
  def N5p(self, p):
    return (1-p[1])*(1+p[2])/8.0
  
  def N5q(self, p):
    return -(1+p[0])*(1+p[2])/8.0
  
  def N5r(self, p):
    return (1+p[0])*(1-p[1])/8.0
  
  def N6(self, p):
    return (1+p[0])*(1+p[1])*(1+p[2])/8.0
  
  def N6p(self, p):
    return (1+p[1])*(1+p[2])/8.0
  
  def N6q(self, p):
    return (1+p[0])*(1+p[2])/8.0
  
  def N6r(self, p):
    return (1+p[0])*(1+p[1])/8.0

  def N7(self, p):
    return (1-p[0])*(1+p[1])*(1+p[2])/8.0
  
  def N7p(self, p):
    return -(1+p[1])*(1+p[2])/8.0
  
  def N7q(self, p):
    return (1-p[0])*(1+p[2])/8.0
  
  def N7r(self, p):
    return (1-p[0])*(1+p[1])/8.0
  
# ----------------------------------------------------------------------
class Hex27(object):

  def __init__(self):
    """
    Setup hex8 cell.
    """
    vertices = numpy.array([[-1.0, -1.0, -1.0], # Corners
                            [+1.0, -1.0, -1.0],
                            [+1.0, +1.0, -1.0],
                            [-1.0, +1.0, -1.0],
                            [-1.0, -1.0, +1.0],
                            [+1.0, -1.0, +1.0],
                            [+1.0, +1.0, +1.0],
                            [-1.0, +1.0, +1.0],
                            [ 0.0, -1.0, -1.0], # Bottom edges
                            [+1.0,  0.0, -1.0],
                            [ 0.0, +1.0, -1.0],
                            [-1.0,  0.0, -1.0],
                            [-1.0, -1.0,  0.0], # Middle edges
                            [+1.0, -1.0,  0.0],
                            [+1.0, +1.0,  0.0],
                            [-1.0, +1.0,  0.0],
                            [ 0.0, -1.0, +1.0], # Top edges
                            [+1.0,  0.0, +1.0],
                            [ 0.0, +1.0, +1.0],
                            [-1.0,  0.0, +1.0],
                            [ 0.0,  0.0,  0.0], # Interior
                            [ 0.0,  0.0, -1.0], # Faces
                            [ 0.0,  0.0, +1.0],
                            [-1.0,  0.0,  0.0],
                            [+1.0,  0.0,  0.0],
                            [ 0.0, -1.0,  0.0],
                            [ 0.0, +1.0,  0.0]])
    quadPts = numpy.array([ [-1.0/3**0.5, -1.0/3**0.5, -1.0/3**0.5],
                            [+1.0/3**0.5, -1.0/3**0.5, -1.0/3**0.5],
                            [+1.0/3**0.5, +1.0/3**0.5, -1.0/3**0.5],
                            [-1.0/3**0.5, +1.0/3**0.5, -1.0/3**0.5],
                            [-1.0/3**0.5, -1.0/3**0.5, +1.0/3**0.5],
                            [+1.0/3**0.5, -1.0/3**0.5, +1.0/3**0.5],
                            [+1.0/3**0.5, +1.0/3**0.5, +1.0/3**0.5],
                            [-1.0/3**0.5, +1.0/3**0.5, +1.0/3**0.5]])
    quadWts = numpy.array( [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    # Compute basis fns and derivatives at quadrature points
    basis = numpy.zeros( (8, 8), dtype=numpy.float64)
    basisDeriv = numpy.zeros( (8, 8, 3), dtype=numpy.float64)
    iQuad = 0
    for q in quadPts:
      basis[iQuad] = numpy.array([self.N0(q), self.N1(q),
                                  self.N2(q), self.N3(q),
                                  self.N4(q), self.N5(q),
                                  self.N6(q), self.N7(q)],
                                 dtype=numpy.float64).reshape( (8,) )
      deriv = numpy.array([[self.N0p(q), self.N0q(q), self.N0r(q)],
                           [self.N1p(q), self.N1q(q), self.N1r(q)],
                           [self.N2p(q), self.N2q(q), self.N2r(q)],
                           [self.N3p(q), self.N3q(q), self.N3r(q)],
                           [self.N4p(q), self.N4q(q), self.N4r(q)],
                           [self.N5p(q), self.N5q(q), self.N5r(q)],
                           [self.N6p(q), self.N6q(q), self.N6r(q)],
                           [self.N7p(q), self.N7q(q), self.N7r(q)]])      
      basisDeriv[iQuad] = deriv.reshape((8, 3))
      iQuad += 1

    self.cellDim = 3
    self.numCorners = len(vertices)
    self.numQuadPts = len(quadPts)
    self.vertices = vertices
    self.quadPts = quadPts
    self.quadWts = quadWts
    self.basis = basis
    self.basisDeriv = basisDeriv
    return


  def N0(self, p):
    return (1-p[0])*(1-p[1])*(1-p[2])/8.0
  
  def N0p(self, p):
    return -(1-p[1])*(1-p[2])/8.0
  
  def N0q(self, p):
    return -(1-p[0])*(1-p[2])/8.0
  
  def N0r(self, p):
    return -(1-p[0])*(1-p[1])/8.0
  
  def N1(self, p):
    return (1+p[0])*(1-p[1])*(1-p[2])/8.0
  
  def N1p(self, p):
    return (1-p[1])*(1-p[2])/8.0
  
  def N1q(self, p):
    return -(1+p[0])*(1-p[2])/8.0
  
  def N1r(self, p):
    return -(1+p[0])*(1-p[1])/8.0
  
  def N2(self, p):
    return (1+p[0])*(1+p[1])*(1-p[2])/8.0
  
  def N2p(self, p):
    return (1+p[1])*(1-p[2])/8.0
  
  def N2q(self, p):
    return (1+p[0])*(1-p[2])/8.0

  def N2r(self, p):
    return -(1+p[0])*(1+p[1])/8.0

  def N3(self, p):
    return (1-p[0])*(1+p[1])*(1-p[2])/8.0
  
  def N3p(self, p):
    return -(1+p[1])*(1-p[2])/8.0
  
  def N3q(self, p):
    return (1-p[0])*(1-p[2])/8.0
  
  def N3r(self, p):
    return -(1-p[0])*(1+p[1])/8.0
  
  def N4(self, p):
    return (1-p[0])*(1-p[1])*(1+p[2])/8.0

  def N4p(self, p):
    return -(1-p[1])*(1+p[2])/8.0
  
  def N4q(self, p):
    return -(1-p[0])*(1+p[2])/8.0
  
  def N4r(self, p):
    return (1-p[0])*(1-p[1])/8.0
  
  def N5(self, p):
    return (1+p[0])*(1-p[1])*(1+p[2])/8.0
  
  def N5p(self, p):
    return (1-p[1])*(1+p[2])/8.0
  
  def N5q(self, p):
    return -(1+p[0])*(1+p[2])/8.0
  
  def N5r(self, p):
    return (1+p[0])*(1-p[1])/8.0
  
  def N6(self, p):
    return (1+p[0])*(1+p[1])*(1+p[2])/8.0
  
  def N6p(self, p):
    return (1+p[1])*(1+p[2])/8.0
  
  def N6q(self, p):
    return (1+p[0])*(1+p[2])/8.0
  
  def N6r(self, p):
    return (1+p[0])*(1+p[1])/8.0

  def N7(self, p):
    return (1-p[0])*(1+p[1])*(1+p[2])/8.0
  
  def N7p(self, p):
    return -(1+p[1])*(1+p[2])/8.0
  
  def N7q(self, p):
    return (1-p[0])*(1+p[2])/8.0
  
  def N7r(self, p):
    return (1-p[0])*(1+p[1])/8.0
  
# ----------------------------------------------------------------------
class TestFIATLagrange(unittest.TestCase):
  """
  Unit testing of FIATLagrange object.
  """

  def test_initialize_line2(self):
    """
    Test initialize() with line2 cell.
    """
    cell = FIATLagrange()
    cell.inventory.dimension = 1
    cell.inventory.degree = 1
    cell.inventory.order  = 1
    cell._configure()
    cell.initialize(spaceDim=1)

    cellE = Line2()
    self._checkVals(cellE, cell)
    from pylith.feassemble.CellGeometry import GeometryLine1D
    self.failUnless(isinstance(cell.geometry, GeometryLine1D))
    return


  def test_initialize_line3(self):
    """
    Test initialize() with line3 cell.
    """
    cell = FIATLagrange()
    cell.inventory.dimension = 1
    cell.inventory.degree = 2
    cell.inventory.order  = 2
    cell._configure()
    cell.initialize(spaceDim=2)

    cellE = Line3()
    self._checkVals(cellE, cell)
    from pylith.feassemble.CellGeometry import GeometryLine2D
    self.failUnless(isinstance(cell.geometry, GeometryLine2D))
    return


  def test_initialize_quad4(self):
    """
    Test initialize() with quad4 cell.
    """
    cell = FIATLagrange()
    cell.inventory.dimension = 2
    cell.inventory.degree = 1
    cell.inventory.order  = 2
    cell._configure()
    cell.initialize(spaceDim=2)

    cellE = Quad4()
    self._checkVals(cellE, cell)
    from pylith.feassemble.CellGeometry import GeometryQuad2D
    self.failUnless(isinstance(cell.geometry, GeometryQuad2D))
    return


  def test_initialize_quad9(self):
    """
    Test initialize() with quad9 cell.
    """
    cell = FIATLagrange()
    cell.inventory.dimension = 2
    cell.inventory.degree = 2
    cell.inventory.order  = 5
    cell._configure()
    cell.initialize(spaceDim=2)

    cellE = Quad9()
    self._checkVals(cellE, cell)
    from pylith.feassemble.CellGeometry import GeometryQuad2D
    self.failUnless(isinstance(cell.geometry, GeometryQuad2D))
    return


  def test_initialize_hex8(self):
    """
    Test initialize() with hex8 cell.
    """
    cell = FIATLagrange()
    cell.inventory.dimension = 3
    cell.inventory.degree = 1
    cell.inventory.order  = 2
    cell._configure()
    cell.initialize(spaceDim=3)

    cellE = Hex8()
    self._checkVals(cellE, cell)
    from pylith.feassemble.CellGeometry import GeometryHex3D
    self.failUnless(isinstance(cell.geometry, GeometryHex3D))
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.feassemble.FIATLagrange import reference_cell
    c = reference_cell()
    return


  def _checkVals(self, cellE, cell):
    """
    Check known values against those generated by FIATLagrange.
    """
    
    # Check basic attributes
    self.assertEqual(cellE.cellDim, cell.cellDim)
    self.assertEqual(cellE.numCorners, cell.numCorners)
    self.assertEqual(cellE.numQuadPts, cell.numQuadPts)

    # Check arrays
    test_double(self, cellE.vertices, cell.vertices)
    test_double(self, cellE.quadPts, cell.quadPts)
    test_double(self, cellE.quadWts, cell.quadWts)
    test_double(self, cellE.basis, cell.basis)
    test_double(self, cellE.basisDeriv, cell.basisDeriv)
    return


# End of file
