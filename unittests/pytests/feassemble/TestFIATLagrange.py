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
# Copyright (c) 2010-2015 University of California, Davis
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
from pylith.utils.testarray import test_scalararray

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
class Quad4Collocated(Quad4):

  def __init__(self):
    """
    Setup quad4 cell.
    """
    vertices = numpy.array([[-1.0, -1.0],
                            [+1.0, -1.0],
                            [+1.0, +1.0],
                            [-1.0, +1.0]])
    quadPts = numpy.array([[-1.0, -1.0],
                           [+1.0, -1.0],
                           [-1.0, +1.0],
                           [+1.0, +1.0]])
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
                            [-1.0, +1.0, -1.0],
                            [+1.0, +1.0, -1.0],
                            [+1.0, -1.0, -1.0],
                            [-1.0, -1.0, +1.0],
                            [+1.0, -1.0, +1.0],
                            [+1.0, +1.0, +1.0],
                            [-1.0, +1.0, +1.0]])
    quadPts = numpy.array([ [-1.0/3**0.5, -1.0/3**0.5, -1.0/3**0.5],
                            [+1.0/3**0.5, -1.0/3**0.5, -1.0/3**0.5],
                            [-1.0/3**0.5, +1.0/3**0.5, -1.0/3**0.5],
                            [+1.0/3**0.5, +1.0/3**0.5, -1.0/3**0.5],
                            [-1.0/3**0.5, -1.0/3**0.5, +1.0/3**0.5],
                            [+1.0/3**0.5, -1.0/3**0.5, +1.0/3**0.5],
                            [-1.0/3**0.5, +1.0/3**0.5, +1.0/3**0.5],
                            [+1.0/3**0.5, +1.0/3**0.5, +1.0/3**0.5]])
    quadWts = numpy.array( [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    # Compute basis fns and derivatives at quadrature points
    basis = numpy.zeros( (8, 8), dtype=numpy.float64)
    basisDeriv = numpy.zeros( (8, 8, 3), dtype=numpy.float64)
    iQuad = 0
    for q in quadPts:
      basis[iQuad] = numpy.array([self.N0(q), self.N3(q),
                                  self.N2(q), self.N1(q),
                                  self.N4(q), self.N5(q),
                                  self.N6(q), self.N7(q)],
                                 dtype=numpy.float64).reshape( (8,) )
      deriv = numpy.array([[self.N0p(q), self.N0q(q), self.N0r(q)],
                           [self.N3p(q), self.N3q(q), self.N3r(q)],
                           [self.N2p(q), self.N2q(q), self.N2r(q)],
                           [self.N1p(q), self.N1q(q), self.N1r(q)],
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
class Hex8Collocated(Hex8):

  def __init__(self):
    """
    Setup hex8 cell.
    """
    vertices = numpy.array([[-1.0, -1.0, -1.0],
                            [-1.0, +1.0, -1.0],
                            [+1.0, +1.0, -1.0],
                            [+1.0, -1.0, -1.0],
                            [-1.0, -1.0, +1.0],
                            [+1.0, -1.0, +1.0],
                            [+1.0, +1.0, +1.0],
                            [-1.0, +1.0, +1.0]])
    quadPts = numpy.array([[-1.0, -1.0, -1.0],
                           [+1.0, -1.0, -1.0],
                           [-1.0, +1.0, -1.0],
                           [+1.0, +1.0, -1.0],
                           [-1.0, -1.0, +1.0],
                           [+1.0, -1.0, +1.0],
                           [-1.0, +1.0, +1.0],
                           [+1.0, +1.0, +1.0]])
    quadWts = numpy.array( [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    # Compute basis fns and derivatives at quadrature points
    basis = numpy.zeros( (8, 8), dtype=numpy.float64)
    basisDeriv = numpy.zeros( (8, 8, 3), dtype=numpy.float64)
    iQuad = 0
    for q in quadPts:
      basis[iQuad] = numpy.array([self.N0(q), self.N3(q),
                                  self.N2(q), self.N1(q),
                                  self.N4(q), self.N5(q),
                                  self.N6(q), self.N7(q)],
                                 dtype=numpy.float64).reshape( (8,) )
      deriv = numpy.array([[self.N0p(q), self.N0q(q), self.N0r(q)],
                           [self.N3p(q), self.N3q(q), self.N3r(q)],
                           [self.N2p(q), self.N2q(q), self.N2r(q)],
                           [self.N1p(q), self.N1q(q), self.N1r(q)],
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


# ----------------------------------------------------------------------
class Hex27(object):

  def __init__(self):
    """
    Setup hex8 cell.
    """
    vertices = numpy.array([[-1.0, -1.0, -1.0], # Corners
                            [-1.0, +1.0, -1.0],
                            [+1.0, +1.0, -1.0],
                            [+1.0, -1.0, -1.0],
                            [-1.0, -1.0, +1.0],
                            [+1.0, -1.0, +1.0],
                            [+1.0, +1.0, +1.0],
                            [-1.0, +1.0, +1.0],
                            [ 0.0, -1.0, -1.0], # Bottom edges
                            [+1.0,  0.0, -1.0],
                            [ 0.0, +1.0, -1.0],
                            [-1.0,  0.0, -1.0],
                            [ 0.0, -1.0, +1.0], # Top edges
                            [+1.0,  0.0, +1.0],
                            [ 0.0, +1.0, +1.0],
                            [-1.0,  0.0, +1.0],
                            [-1.0, -1.0,  0.0], # Middle edges
                            [+1.0, -1.0,  0.0],
                            [+1.0, +1.0,  0.0],
                            [-1.0, +1.0,  0.0],
                            [-1.0,  0.0,  0.0], # Faces
                            [+1.0,  0.0,  0.0],
                            [ 0.0, -1.0,  0.0],
                            [ 0.0, +1.0,  0.0],
                            [ 0.0,  0.0, -1.0],
                            [ 0.0,  0.0, +1.0],
                            [ 0.0,  0.0,  0.0]]) # Interior
    quadPts = numpy.array([ [-(3.0/5)**0.5, -(3.0/5)**0.5, -(3.0/5)**0.5],
                            [          0.0, -(3.0/5)**0.5, -(3.0/5)**0.5],
                            [+(3.0/5)**0.5, -(3.0/5)**0.5, -(3.0/5)**0.5],
                            [-(3.0/5)**0.5,           0.0, -(3.0/5)**0.5],
                            [          0.0,           0.0, -(3.0/5)**0.5],
                            [+(3.0/5)**0.5,           0.0, -(3.0/5)**0.5],
                            [-(3.0/5)**0.5, +(3.0/5)**0.5, -(3.0/5)**0.5],
                            [          0.0, +(3.0/5)**0.5, -(3.0/5)**0.5],
                            [+(3.0/5)**0.5, +(3.0/5)**0.5, -(3.0/5)**0.5],
                            [-(3.0/5)**0.5, -(3.0/5)**0.5,           0.0],
                            [          0.0, -(3.0/5)**0.5,           0.0],
                            [+(3.0/5)**0.5, -(3.0/5)**0.5,           0.0],
                            [-(3.0/5)**0.5,           0.0,           0.0],
                            [          0.0,           0.0,           0.0],
                            [+(3.0/5)**0.5,           0.0,           0.0],
                            [-(3.0/5)**0.5, +(3.0/5)**0.5,           0.0],
                            [          0.0, +(3.0/5)**0.5,           0.0],
                            [+(3.0/5)**0.5, +(3.0/5)**0.5,           0.0],
                            [-(3.0/5)**0.5, -(3.0/5)**0.5, +(3.0/5)**0.5],
                            [          0.0, -(3.0/5)**0.5, +(3.0/5)**0.5],
                            [+(3.0/5)**0.5, -(3.0/5)**0.5, +(3.0/5)**0.5],
                            [-(3.0/5)**0.5,           0.0, +(3.0/5)**0.5],
                            [          0.0,           0.0, +(3.0/5)**0.5],
                            [+(3.0/5)**0.5,           0.0, +(3.0/5)**0.5],
                            [-(3.0/5)**0.5, +(3.0/5)**0.5, +(3.0/5)**0.5],
                            [          0.0, +(3.0/5)**0.5, +(3.0/5)**0.5],
                            [+(3.0/5)**0.5, +(3.0/5)**0.5, +(3.0/5)**0.5] ])
    quadWts = numpy.array( [125.0/729, 200.0/729, 125.0/729,
                            200.0/729, 320.0/729, 200.0/729,
                            125.0/729, 200.0/729, 125.0/729,
                            200.0/729, 320.0/729, 200.0/729,
                            320.0/729, 512.0/729, 320.0/729,
                            200.0/729, 320.0/729, 200.0/729,
                            125.0/729, 200.0/729, 125.0/729,
                            200.0/729, 320.0/729, 200.0/729,
                            125.0/729, 200.0/729, 125.0/729])

    # Compute basis fns and derivatives at quadrature points
    basis = numpy.zeros( (27, 27), dtype=numpy.float64)
    basisDeriv = numpy.zeros( (27, 27, 3), dtype=numpy.float64)
    iQuad = 0
    for q in quadPts:
      basis[iQuad] = numpy.array([self.N0(q), self.N3(q), self.N2(q), self.N1(q), # Corners
                                  self.N4(q), self.N5(q), self.N6(q), self.N7(q),
                                  self.N8(q), self.N9(q), self.N10(q), self.N11(q), # Edges
                                  self.N12(q), self.N13(q), self.N14(q), self.N15(q),
                                  self.N16(q), self.N17(q), self.N18(q), self.N19(q),
                                  self.N20(q), # Interior
                                  self.N21(q), self.N22(q), # Faces
                                  self.N23(q), self.N24(q),
                                  self.N25(q), self.N26(q)],
                                 dtype=numpy.float64).reshape( (27,) )
      deriv = numpy.array([[self.N0p(q), self.N0q(q), self.N0r(q)],
                           [self.N3p(q), self.N3q(q), self.N3r(q)],
                           [self.N2p(q), self.N2q(q), self.N2r(q)],
                           [self.N1p(q), self.N1q(q), self.N1r(q)],
                           [self.N4p(q), self.N4q(q), self.N4r(q)],
                           [self.N5p(q), self.N5q(q), self.N5r(q)],
                           [self.N6p(q), self.N6q(q), self.N6r(q)],
                           [self.N7p(q), self.N7q(q), self.N7r(q)],
                           [self.N8p(q), self.N8q(q), self.N8r(q)],
                           [self.N9p(q), self.N9q(q), self.N9r(q)],
                           [self.N10p(q), self.N10q(q), self.N10r(q)],
                           [self.N11p(q), self.N11q(q), self.N11r(q)],
                           [self.N12p(q), self.N12q(q), self.N12r(q)],
                           [self.N13p(q), self.N13q(q), self.N13r(q)],
                           [self.N14p(q), self.N14q(q), self.N14r(q)],
                           [self.N15p(q), self.N15q(q), self.N15r(q)],
                           [self.N16p(q), self.N16q(q), self.N16r(q)],
                           [self.N17p(q), self.N17q(q), self.N17r(q)],
                           [self.N18p(q), self.N18q(q), self.N18r(q)],
                           [self.N19p(q), self.N19q(q), self.N19r(q)],
                           [self.N20p(q), self.N20q(q), self.N20r(q)],
                           [self.N21p(q), self.N21q(q), self.N21r(q)],
                           [self.N22p(q), self.N22q(q), self.N22r(q)],
                           [self.N23p(q), self.N23q(q), self.N23r(q)],
                           [self.N24p(q), self.N24q(q), self.N24r(q)],
                           [self.N25p(q), self.N25q(q), self.N25r(q)],
                           [self.N26p(q), self.N26q(q), self.N26r(q)]])      
      basisDeriv[iQuad] = deriv.reshape((27, 3))
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


  # Corners
  def N0(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(-0.5*p[1])*(1.0-p[1])*(-0.5*p[2])*(1.0-p[2])

  def N0p(self, p):
    return (p[0]-0.5)*(-0.5*p[1])*(1.0-p[1])*(-0.5*p[2])*(1.0-p[2])

  def N0q(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(p[1]-0.5)*(-0.5*p[2])*(1.0-p[2])

  def N0r(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(-0.5*p[1])*(1.0-p[1])*(p[2]-0.5)


  def N1(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(-0.5*p[1])*(1.0-p[1])*(-0.5*p[2])*(1.0-p[2])

  def N1p(self, p):
    return (p[0]+0.5)*(-0.5*p[1])*(1.0-p[1])*(-0.5*p[2])*(1.0-p[2])

  def N1q(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(p[1]-0.5)*(-0.5*p[2])*(1.0-p[2])

  def N1r(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(-0.5*p[1])*(1.0-p[1])*(p[2]-0.5)


  def N2(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(+0.5*p[1])*(1.0+p[1])*(-0.5*p[2])*(1.0-p[2])

  def N2p(self, p):
    return (p[0]+0.5)*(+0.5*p[1])*(1.0+p[1])*(-0.5*p[2])*(1.0-p[2])

  def N2q(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(p[1]+0.5)*(-0.5*p[2])*(1.0-p[2])

  def N2r(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(+0.5*p[1])*(1.0+p[1])*(p[2]-0.5)


  def N3(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(+0.5*p[1])*(1.0+p[1])*(-0.5*p[2])*(1.0-p[2])

  def N3p(self, p):
    return (p[0]-0.5)*(+0.5*p[1])*(1.0+p[1])*(-0.5*p[2])*(1.0-p[2])

  def N3q(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(p[1]+0.5)*(-0.5*p[2])*(1.0-p[2])

  def N3r(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(+0.5*p[1])*(1.0+p[1])*(p[2]-0.5)


  def N4(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(-0.5*p[1])*(1.0-p[1])*(+0.5*p[2])*(1.0+p[2])

  def N4p(self, p):
    return (p[0]-0.5)*(-0.5*p[1])*(1.0-p[1])*(+0.5*p[2])*(1.0+p[2])

  def N4q(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(p[1]-0.5)*(+0.5*p[2])*(1.0+p[2])

  def N4r(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(-0.5*p[1])*(1.0-p[1])*(p[2]+0.5)


  def N5(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(-0.5*p[1])*(1.0-p[1])*(+0.5*p[2])*(1.0+p[2])

  def N5p(self, p):
    return (p[0]+0.5)*(-0.5*p[1])*(1.0-p[1])*(+0.5*p[2])*(1.0+p[2])

  def N5q(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(p[1]-0.5)*(+0.5*p[2])*(1.0+p[2])

  def N5r(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(-0.5*p[1])*(1.0-p[1])*(p[2]+0.5)


  def N6(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(+0.5*p[1])*(1.0+p[1])*(+0.5*p[2])*(1.0+p[2])

  def N6p(self, p):
    return (p[0]+0.5)*(+0.5*p[1])*(1.0+p[1])*(+0.5*p[2])*(1.0+p[2])

  def N6q(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(p[1]+0.5)*(+0.5*p[2])*(1.0+p[2])

  def N6r(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(+0.5*p[1])*(1.0+p[1])*(p[2]+0.5)


  def N7(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(+0.5*p[1])*(1.0+p[1])*(+0.5*p[2])*(1.0+p[2])

  def N7p(self, p):
    return (p[0]-0.5)*(+0.5*p[1])*(1.0+p[1])*(+0.5*p[2])*(1.0+p[2])

  def N7q(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(p[1]+0.5)*(+0.5*p[2])*(1.0+p[2])

  def N7r(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(+0.5*p[1])*(1.0+p[1])*(p[2]+0.5)


  # Bottom edges
  def N8(self, p):
    return (1.0-p[0]**2)*(-0.5*p[1])*(1.0-p[1])*(-0.5*p[2])*(1.0-p[2])

  def N8p(self, p):
    return (-2.0*p[0])*(-0.5*p[1])*(1.0-p[1])*(-0.5*p[2])*(1.0-p[2])

  def N8q(self, p):
    return (1.0-p[0]**2)*(p[1]-0.5)*(-0.5*p[2])*(1.0-p[2])

  def N8r(self, p):
    return (1.0-p[0]**2)*(-0.5*p[1])*(1.0-p[1])*(p[2]-0.5)


  def N9(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(1.0-p[1]**2)*(-0.5*p[2])*(1.0-p[2])

  def N9p(self, p):
    return (p[0]+0.5)*(1.0-p[1]**2)*(-0.5*p[2])*(1.0-p[2])

  def N9q(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(-2.0*p[1])*(-0.5*p[2])*(1.0-p[2])

  def N9r(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(1.0-p[1]**2)*(p[2]-0.5)


  def N10(self, p):
    return (1.0-p[0]**2)*(+0.5*p[1])*(1.0+p[1])*(-0.5*p[2])*(1.0-p[2])

  def N10p(self, p):
    return (-2.0*p[0])*(+0.5*p[1])*(1.0+p[1])*(-0.5*p[2])*(1.0-p[2])

  def N10q(self, p):
    return (1.0-p[0]**2)*(p[1]+0.5)*(-0.5*p[2])*(1.0-p[2])

  def N10r(self, p):
    return (1.0-p[0]**2)*(+0.5*p[1])*(1.0+p[1])*(p[2]-0.5)


  def N11(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(1.0-p[1]**2)*(-0.5*p[2])*(1.0-p[2])

  def N11p(self, p):
    return (p[0]-0.5)*(1.0-p[1]**2)*(-0.5*p[2])*(1.0-p[2])

  def N11q(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(-2.0*p[1])*(-0.5*p[2])*(1.0-p[2])

  def N11r(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(1.0-p[1]**2)*(p[2]-0.5)


  # Top edges
  def N12(self, p):
    return (1.0-p[0]**2)*(-0.5*p[1])*(1.0-p[1])*(+0.5*p[2])*(1.0+p[2])

  def N12p(self, p):
    return (-2.0*p[0])*(-0.5*p[1])*(1.0-p[1])*(+0.5*p[2])*(1.0+p[2])

  def N12q(self, p):
    return (1.0-p[0]**2)*(p[1]-0.5)*(+0.5*p[2])*(1.0+p[2])

  def N12r(self, p):
    return (1.0-p[0]**2)*(-0.5*p[1])*(1.0-p[1])*(p[2]+0.5)


  def N13(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(1.0-p[1]**2)*(+0.5*p[2])*(1.0+p[2])

  def N13p(self, p):
    return (p[0]+0.5)*(1.0-p[1]**2)*(+0.5*p[2])*(1.0+p[2])

  def N13q(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(-2.0*p[1])*(+0.5*p[2])*(1.0+p[2])

  def N13r(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(1.0-p[1]**2)*(p[2]+0.5)


  def N14(self, p):
    return (1.0-p[0]**2)*(+0.5*p[1])*(1.0+p[1])*(+0.5*p[2])*(1.0+p[2])

  def N14p(self, p):
    return (-2.0*p[0])*(+0.5*p[1])*(1.0+p[1])*(+0.5*p[2])*(1.0+p[2])

  def N14q(self, p):
    return (1.0-p[0]**2)*(p[1]+0.5)*(+0.5*p[2])*(1.0+p[2])

  def N14r(self, p):
    return (1.0-p[0]**2)*(+0.5*p[1])*(1.0+p[1])*(p[2]+0.5)


  def N15(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(1.0-p[1]**2)*(+0.5*p[2])*(1.0+p[2])

  def N15p(self, p):
    return (p[0]-0.5)*(1.0-p[1]**2)*(+0.5*p[2])*(1.0+p[2])

  def N15q(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(-2.0*p[1])*(+0.5*p[2])*(1.0+p[2])

  def N15r(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(1.0-p[1]**2)*(p[2]+0.5)


  # Middle edges
  def N16(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(-0.5*p[1])*(1.0-p[1])*(1.0-p[2]**2)

  def N16p(self, p):
    return (p[0]-0.5)*(-0.5*p[1])*(1.0-p[1])*(1.0-p[2]**2)

  def N16q(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(p[1]-0.5)*(1.0-p[2]**2)

  def N16r(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(-0.5*p[1])*(1.0-p[1])*(-2.0*p[2])


  def N17(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(-0.5*p[1])*(1.0-p[1])*(1.0-p[2]**2)

  def N17p(self, p):
    return (p[0]+0.5)*(-0.5*p[1])*(1.0-p[1])*(1.0-p[2]**2)

  def N17q(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(p[1]-0.5)*(1.0-p[2]**2)

  def N17r(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(-0.5*p[1])*(1.0-p[1])*(-2.0*p[2])


  def N18(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(+0.5*p[1])*(1.0+p[1])*(1.0-p[2]**2)

  def N18p(self, p):
    return (p[0]+0.5)*(+0.5*p[1])*(1.0+p[1])*(1.0-p[2]**2)

  def N18q(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(p[1]+0.5)*(1.0-p[2]**2)

  def N18r(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(+0.5*p[1])*(1.0+p[1])*(-2.0*p[2])


  def N19(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(+0.5*p[1])*(1.0+p[1])*(1.0-p[2]**2)

  def N19p(self, p):
    return (p[0]-0.5)*(+0.5*p[1])*(1.0+p[1])*(1.0-p[2]**2)

  def N19q(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(p[1]+0.5)*(1.0-p[2]**2)

  def N19r(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(+0.5*p[1])*(1.0+p[1])*(-2.0*p[2])


  # Left/right
  def N20(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(1.0-p[1]**2)*(1.0-p[2]**2)

  def N20p(self, p):
    return (p[0]-0.5)*(1.0-p[1]**2)*(1.0-p[2]**2)

  def N20q(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(-2.0*p[1])*(1.0-p[2]**2)

  def N20r(self, p):
    return (-0.5*p[0])*(1.0-p[0])*(1.0-p[1]**2)*(-2.0*p[2])


  def N21(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(1.0-p[1]**2)*(1.0-p[2]**2)

  def N21p(self, p):
    return (p[0]+0.5)*(1.0-p[1]**2)*(1.0-p[2]**2)

  def N21q(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(-2.0*p[1])*(1.0-p[2]**2)

  def N21r(self, p):
    return (+0.5*p[0])*(1.0+p[0])*(1.0-p[1]**2)*(-2.0*p[2])


  # Front/back
  def N22(self, p):
    return (1.0-p[0]**2)*(-0.5*p[1])*(1.0-p[1])*(1.0-p[2]**2)

  def N22p(self, p):
    return (-2.0*p[0])*(-0.5*p[1])*(1.0-p[1])*(1.0-p[2]**2)

  def N22q(self, p):
    return (1.0-p[0]**2)*(p[1]-0.5)*(1.0-p[2]**2)

  def N22r(self, p):
    return (1.0-p[0]**2)*(-0.5*p[1])*(1.0-p[1])*(-2.0*p[2])


  def N23(self, p):
    return (1.0-p[0]**2)*(+0.5*p[1])*(1.0+p[1])*(1.0-p[2]**2)

  def N23p(self, p):
    return (-2.0*p[0])*(+0.5*p[1])*(1.0+p[1])*(1.0-p[2]**2)

  def N23q(self, p):
    return (1.0-p[0]**2)*(p[1]+0.5)*(1.0-p[2]**2)

  def N23r(self, p):
    return (1.0-p[0]**2)*(+0.5*p[1])*(1.0+p[1])*(-2.0*p[2])


  # Bottom/top
  def N24(self, p):
    return (1.0-p[0]**2)*(1.0-p[1]**2)*(-0.5*p[2])*(1.0-p[2])

  def N24p(self, p):
    return (-2.0*p[0])*(1.0-p[1]**2)*(-0.5*p[2])*(1.0-p[2])

  def N24q(self, p):
    return (1.0-p[0]**2)*(-2.0*p[1])*(-0.5*p[2])*(1.0-p[2])

  def N24r(self, p):
    return (1.0-p[0]**2)*(1.0-p[1]**2)*(p[2]-0.5)


  def N25(self, p):
    return (1.0-p[0]**2)*(1.0-p[1]**2)*(+0.5*p[2])*(1.0+p[2])

  def N25p(self, p):
    return (-2.0*p[0])*(1.0-p[1]**2)*(+0.5*p[2])*(1.0+p[2])

  def N25q(self, p):
    return (1.0-p[0]**2)*(-2.0*p[1])*(+0.5*p[2])*(1.0+p[2])

  def N25r(self, p):
    return (1.0-p[0]**2)*(1.0-p[1]**2)*(p[2]+0.5)


  # Interior
  def N26(self, p):
    return (1.0-p[0]**2)*(1.0-p[1]**2)*(1.0-p[2]**2)

  def N26p(self, p):
    return (-2.0*p[0])*(1.0-p[1]**2)*(1.0-p[2]**2)

  def N26q(self, p):
    return (1.0-p[0]**2)*(-2.0*p[1])*(1.0-p[2]**2)

  def N26r(self, p):
    return (1.0-p[0]**2)*(1.0-p[1]**2)*(-2.0*p[2])


# ----------------------------------------------------------------------
class TestFIATLagrange(unittest.TestCase):
  """
  Unit testing of FIATLagrange object.
  """

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


  def test_initialize_quad4_collocated(self):
    """
    Test initialize() with quad4 cell.
    """
    cell = FIATLagrange()
    cell.inventory.dimension = 2
    cell.inventory.degree = 1
    cell.inventory.order  = 1
    cell.inventory.collocateQuad = True
    cell._configure()
    cell.initialize(spaceDim=2)

    cellE = Quad4Collocated()
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
    cell.inventory.order  = 3
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


  def test_initialize_hex8_collocated(self):
    """
    Test initialize() with hex8 cell.
    """
    cell = FIATLagrange()
    cell.inventory.dimension = 3
    cell.inventory.degree = 1
    cell.inventory.order  = 1
    cell.inventory.collocateQuad = True
    cell._configure()
    cell.initialize(spaceDim=3)

    cellE = Hex8Collocated()
    self._checkVals(cellE, cell)
    from pylith.feassemble.CellGeometry import GeometryHex3D
    self.failUnless(isinstance(cell.geometry, GeometryHex3D))
    return


  def test_initialize_hex27(self):
    """
    Test initialize() with hex27 cell.
    """
    cell = FIATLagrange()
    cell.inventory.dimension = 3
    cell.inventory.degree = 2
    cell.inventory.order  = 3
    cell._configure()
    cell.initialize(spaceDim=3)

    cellE = Hex27()
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
    test_scalararray(self, cellE.vertices, cell.vertices)
    test_scalararray(self, cellE.quadPts, cell.quadPts)
    test_scalararray(self, cellE.quadWts, cell.quadWts)
    test_scalararray(self, cellE.basis, cell.basis)
    test_scalararray(self, cellE.basisDeriv, cell.basisDeriv)
    return


# End of file
