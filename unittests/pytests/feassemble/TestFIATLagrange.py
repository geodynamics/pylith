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
                            [-1.0, +1.0],
                            [+1.0, +1.0]])
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
    return (1-p[0])*(1+p[1])/4.0

  def N2p(self, p):
    return -(1+p[1])/4.0

  def N2q(self, p):
    return (1-p[0])/4.0

  def N3(self, p):
    return (1+p[0])*(1+p[1])/4.0

  def N3p(self, p):
    return (1+p[1])/4.0

  def N3q(self, p):
    return (1+p[0])/4.0

# ----------------------------------------------------------------------
class Hex8(object):

  def __init__(self):
    """
    Setup hex8 cell.
    """
    vertices = numpy.array([[-1.0, -1.0, -1.0],
                            [+1.0, -1.0, -1.0],
                            [-1.0, +1.0, -1.0],
                            [+1.0, +1.0, -1.0],
                            [-1.0, -1.0, +1.0],
                            [+1.0, -1.0, +1.0],
                            [-1.0, +1.0, +1.0],
                            [+1.0, +1.0, +1.0]])
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
    return (1-p[0])*(1+p[1])*(1-p[2])/8.0
  
  def N2p(self, p):
    return -(1+p[1])*(1-p[2])/8.0
  
  def N2q(self, p):
    return (1-p[0])*(1-p[2])/8.0
  
  def N2r(self, p):
    return -(1-p[0])*(1+p[1])/8.0
  
  def N3(self, p):
    return (1+p[0])*(1+p[1])*(1-p[2])/8.0
  
  def N3p(self, p):
    return (1+p[1])*(1-p[2])/8.0
  
  def N3q(self, p):
    return (1+p[0])*(1-p[2])/8.0

  def N3r(self, p):
    return -(1+p[0])*(1+p[1])/8.0

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
    return (1-p[0])*(1+p[1])*(1+p[2])/8.0
  
  def N6p(self, p):
    return -(1+p[1])*(1+p[2])/8.0
  
  def N6q(self, p):
    return (1-p[0])*(1+p[2])/8.0
  
  def N6r(self, p):
    return (1-p[0])*(1+p[1])/8.0
  
  def N7(self, p):
    return (1+p[0])*(1+p[1])*(1+p[2])/8.0
  
  def N7p(self, p):
    return (1+p[1])*(1+p[2])/8.0
  
  def N7q(self, p):
    return (1+p[0])*(1+p[2])/8.0
  
  def N7r(self, p):
    return (1+p[0])*(1+p[1])/8.0

# ----------------------------------------------------------------------
class TestFIATLagrange(unittest.TestCase):
  """
  Unit testing of FIATLagrange object.
  """

  def test_initialize_line2(self):
    """
    Test initialize() with line2 cell.
    """
    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATLagrange()
    cell.cellDim = 1
    cell.degree = 1
    cell.order  = 1
    cell.initialize(spaceDim=1)

    cellE = Line2()
    self._checkVals(cellE, cell)
    from pylith.feassemble.geometry.GeometryLine1D import GeometryLine1D
    self.assertTrue(isinstance(cell.geometry, GeometryLine1D))
    return


  def test_initialize_line3(self):
    """
    Test initialize() with line3 cell.
    """
    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATLagrange()
    cell.cellDim = 1
    cell.degree = 2
    cell.order  = 2
    cell.initialize(spaceDim=2)

    cellE = Line3()
    self._checkVals(cellE, cell)
    from pylith.feassemble.geometry.GeometryLine2D import GeometryLine2D
    self.assertTrue(isinstance(cell.geometry, GeometryLine2D))
    return


  def test_initialize_quad4(self):
    """
    Test initialize() with quad4 cell.
    """
    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATLagrange()
    cell.cellDim = 2
    cell.degree = 1
    cell.order  = 2
    cell.initialize(spaceDim=2)

    cellE = Quad4()
    self._checkVals(cellE, cell)
    from pylith.feassemble.geometry.GeometryQuad2D import GeometryQuad2D
    self.assertTrue(isinstance(cell.geometry, GeometryQuad2D))
    return


  def test_initialize_hex8(self):
    """
    Test initialize() with hex8 cell.
    """
    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATLagrange()
    cell.cellDim = 3
    cell.degree = 1
    cell.order  = 2
    cell.initialize(spaceDim=3)

    cellE = Hex8()
    self._checkVals(cellE, cell)
    from pylith.feassemble.geometry.GeometryHex3D import GeometryHex3D
    self.assertTrue(isinstance(cell.geometry, GeometryHex3D))
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
