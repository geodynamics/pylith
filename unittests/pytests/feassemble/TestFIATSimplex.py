#!/usr/bin/env python
#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# {LicenseText}
#
# ======================================================================
#

## @file unittests/pytests/feassemble/TestFIATSimplex.py

## @brief Unit testing of FIATSimplex object.

import unittest
import numpy
from pylith.feassemble.FIATSimplex import FIATSimplex
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
class Tri3(object):

  def __init__(self):
    """
    Setup tri33 cell.
    """
    vertices = numpy.array([[-1.0, -1.0],
                            [+1.0, -1.0],
                            [-1.0, +1.0]])
    quadPts = numpy.array([ [-1.0/3.0, -1.0/3.0] ])
    quadWts = numpy.array( [2.0])

    # Compute basis fns and derivatives at quadrature points
    basis = numpy.zeros( (1, 3), dtype=numpy.float64)
    basisDeriv = numpy.zeros( (1, 3, 2), dtype=numpy.float64)
    iQuad = 0
    for q in quadPts:
      basis[iQuad] = numpy.array([self.N0(q), self.N1(q), self.N2(q)],
                                 dtype=numpy.float64).reshape( (3,) )
      deriv = numpy.array([[self.N0p(q), self.N0q(q)],
                           [self.N1p(q), self.N1q(q)],
                           [self.N2p(q), self.N2q(q)]])      
      basisDeriv[iQuad] = deriv.reshape((3, 2))
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
    return 0.5*(-p[0]-p[1])

  def N0p(self, p):
    return -0.5

  def N0q(self, p):
    return -0.5

  def N1(self, p):
    return 0.5*(1.0+p[0])

  def N1p(self, p):
    return 0.5

  def N1q(self, p):
    return 0.0

  def N2(self, p):
    return 0.5*(1.0+p[1])

  def N2p(self, p):
    return 0.0

  def N2q(self, p):
    return 0.5

# ----------------------------------------------------------------------
class Tet4(object):

  def __init__(self):
    """
    Setup tri33 cell.
    """
    vertices = numpy.array([[-1.0, -1.0, -1.0],
                            [+1.0, -1.0, -1.0],
                            [-1.0, +1.0, -1.0],
                            [-1.0, -1.0, +1.0]])
    quadPts = numpy.array([ [-1.0/2.0, -1.0/2.0, -1.0/2.0] ])
    quadWts = numpy.array( [4.0/3.0])

    # Compute basis fns and derivatives at quadrature points
    basis = numpy.zeros( (1, 4), dtype=numpy.float64)
    basisDeriv = numpy.zeros( (1, 4, 3), dtype=numpy.float64)
    iQuad = 0
    for q in quadPts:
      basis[iQuad] = numpy.array([self.N0(q), self.N1(q), 
                                  self.N2(q), self.N3(q)],
                                 dtype=numpy.float64).reshape( (4,) )
      deriv = numpy.array([[self.N0p(q), self.N0q(q), self.N0r(q)],
                           [self.N1p(q), self.N1q(q), self.N1r(q)],
                           [self.N2p(q), self.N2q(q), self.N2r(q)],
                           [self.N3p(q), self.N3q(q), self.N3r(q)]])      
      basisDeriv[iQuad] = deriv.reshape((4, 3))
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
    return 0.5*(-1-p[0]-p[1]-p[2])

  def N0p(self, p):
    return -0.5

  def N0q(self, p):
    return -0.5

  def N0r(self, p):
    return -0.5

  def N1(self, p):
    return 0.5*(1.0+p[0])

  def N1p(self, p):
    return 0.5

  def N1q(self, p):
    return 0.0

  def N1r(self, p):
    return 0.0

  def N2(self, p):
    return 0.5*(1.0+p[1])

  def N2p(self, p):
    return 0.0

  def N2q(self, p):
    return 0.5

  def N2r(self, p):
    return 0.0

  def N3(self, p):
    return 0.5*(1.0+p[2])

  def N3p(self, p):
    return 0.0

  def N3q(self, p):
    return 0.0

  def N3r(self, p):
    return 0.5

# ----------------------------------------------------------------------
class TestFIATSimplex(unittest.TestCase):
  """
  Unit testing of FIATSimplex object.
  """
  

  def test_shape(self):
    """
    Test _getShape().
    """
    cell = FIATSimplex()

    import FIAT.shapes

    cell.shape = "line"
    shape = cell._getShape()
    self.assertEqual(FIAT.shapes.LINE, shape)

    cell.shape = "triangle"
    shape = cell._getShape()
    self.assertEqual(FIAT.shapes.TRIANGLE, shape)

    cell.shape = "tetrahedron"
    shape = cell._getShape()
    self.assertEqual(FIAT.shapes.TETRAHEDRON, shape)
    return


  def test_initialize_line2(self):
    """
    Test initialize() with line2 cell.
    """
    cell = FIATSimplex()
    cell.shape  = "line"
    cell.degree = 1
    cell.order  = 1
    cell.initialize()

    cellE = Line2()
    self._checkVals(cellE, cell)
    return


  def test_initialize_line3(self):
    """
    Test initialize() with line3 cell.
    """
    cell = FIATSimplex()
    cell.shape  = "line"
    cell.degree = 2
    cell.order  = 2
    cell.initialize()

    cellE = Line3()
    self._checkVals(cellE, cell)
    return


  def test_initialize_tri3(self):
    """
    Test initialize() with tri3 cell.
    """
    cell = FIATSimplex()
    cell.shape  = "triangle"
    cell.degree = 1
    cell.order  = 1
    cell.initialize()

    cellE = Tri3()
    self._checkVals(cellE, cell)
    return


  def test_initialize_tet4(self):
    """
    Test initialize() with tet4 cell.
    """
    cell = FIATSimplex()
    cell.shape  = "tetrahedron"
    cell.degree = 1
    cell.order  = 1
    cell.initialize()

    cellE = Tet4()
    self._checkVals(cellE, cell)
    return


  def _checkVals(self, cellE, cell):
    """
    Check known values against those generated by FIATSimplex.
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
