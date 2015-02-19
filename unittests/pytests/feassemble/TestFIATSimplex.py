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

## @file unittests/pytests/feassemble/TestFIATSimplex.py

## @brief Unit testing of FIATSimplex object.

import unittest
import numpy
from pylith.feassemble.FIATSimplex import FIATSimplex
from pylith.utils.testarray import test_scalararray

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
class Tri3Collocated(Tri3):

  def __init__(self):
    """
    Setup tri33 cell.
    """
    vertices = numpy.array([[-1.0, -1.0],
                            [+1.0, -1.0],
                            [-1.0, +1.0]])
    quadPts = vertices[:]
    quadWts = numpy.array( [2.0/3.0, 2.0/3.0, 2.0/3.0])

    # Compute basis fns and derivatives at quadrature points
    basis = numpy.zeros( (3, 3), dtype=numpy.float64)
    basisDeriv = numpy.zeros( (3, 3, 2), dtype=numpy.float64)
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


# ----------------------------------------------------------------------
class Tri6(object):

  def __init__(self):
    """
    Setup tri33 cell.
    """
    vertices = numpy.array([[-1.0, -1.0],
                            [+1.0, -1.0],
                            [-1.0, +1.0],
                            [ 0.0, -1.0],
                            [ 0.0,  0.0],
                            [-1.0,  0.0]])
    quadPts = numpy.array([ [-0.64288254, -0.68989795],
                            [-0.84993778,  0.28989795],
                            [ 0.33278049, -0.68989795],
                            [-0.43996017,  0.28989795]])
    quadWts = numpy.array( [0.63608276,  0.36391724,  0.63608276,  0.36391724])

    # Compute basis fns and derivatives at quadrature points
    basis = numpy.zeros( (4, 6), dtype=numpy.float64)
    basisDeriv = numpy.zeros( (4, 6, 2), dtype=numpy.float64)
    iQuad = 0
    for q in quadPts:
      basis[iQuad] = numpy.array([self.N0(q), self.N1(q), self.N2(q),
                                  self.N3(q), self.N4(q), self.N5(q)],
                                 dtype=numpy.float64).reshape( (6,) )
      deriv = numpy.array([[self.N0p(q), self.N0q(q)],
                           [self.N1p(q), self.N1q(q)],
                           [self.N2p(q), self.N2q(q)],
                           [self.N3p(q), self.N3q(q)],
                           [self.N4p(q), self.N4q(q)],
                           [self.N5p(q), self.N5q(q)]])
      basisDeriv[iQuad] = deriv.reshape((6, 2))
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
    return 0.5*(-p[0]-p[1])*(-1.0-p[0]-p[1])

  def N0p(self, p):
    return 0.5+p[0]+p[1]

  def N0q(self, p):
    return 0.5+p[0]+p[1]


  def N1(self, p):
    return 0.5*(1.0+p[0])*(p[0])

  def N1p(self, p):
    return 0.5+p[0]

  def N1q(self, p):
    return 0


  def N2(self, p):
    return 0.5*(1.0+p[1])*(p[1])

  def N2p(self, p):
    return 0

  def N2q(self, p):
    return 0.5+p[1]


  def N3(self, p):
    return (-p[0]-p[1])*(1+p[0])

  def N3p(self, p):
    return -1.0-2*p[0]-p[1]

  def N3q(self, p):
    return -(1+p[0])


  def N4(self, p):
    return (1.0+p[0])*(1+p[1])

  def N4p(self, p):
    return (1+p[1])

  def N4q(self, p):
    return (1.0+p[0])


  def N5(self, p):
    return (-p[0]-p[1])*(1+p[1])

  def N5p(self, p):
    return -(1+p[1])

  def N5q(self, p):
    return -1.0-p[0]-2*p[1]


# ----------------------------------------------------------------------
class Tet4(object):

  def __init__(self):
    """
    Setup tri33 cell.
    """
    vertices = numpy.array([[+1.0, -1.0, -1.0],
                            [-1.0, -1.0, -1.0],
                            [-1.0, +1.0, -1.0],
                            [-1.0, -1.0, +1.0]])
    quadPts = numpy.array([ [-1.0/2.0, -1.0/2.0, -1.0/2.0] ])
    quadWts = numpy.array( [4.0/3.0])

    # Compute basis fns and derivatives at quadrature points
    basis = numpy.zeros( (1, 4), dtype=numpy.float64)
    basisDeriv = numpy.zeros( (1, 4, 3), dtype=numpy.float64)
    iQuad = 0
    for q in quadPts:
      basis[iQuad] = numpy.array([self.N1(q), self.N0(q), 
                                  self.N2(q), self.N3(q)],
                                 dtype=numpy.float64).reshape( (4,) )
      deriv = numpy.array([[self.N1p(q), self.N1q(q), self.N1r(q)],
                           [self.N0p(q), self.N0q(q), self.N0r(q)],
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
class Tet4Collocated(Tet4):

  def __init__(self):
    """
    Setup tri33 cell.
    """
    vertices = numpy.array([[+1.0, -1.0, -1.0],
                            [-1.0, -1.0, -1.0],
                            [-1.0, +1.0, -1.0],
                            [-1.0, -1.0, +1.0]])
    quadPts = numpy.array([[-1.0, -1.0, -1.0],
                            [+1.0, -1.0, -1.0],
                            [-1.0, +1.0, -1.0],
                            [-1.0, -1.0, +1.0]])
    quadWts = numpy.array( [1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/3.0])

    # Compute basis fns and derivatives at quadrature points
    basis = numpy.zeros( (4, 4), dtype=numpy.float64)
    basisDeriv = numpy.zeros( (4, 4, 3), dtype=numpy.float64)
    iQuad = 0
    for q in quadPts:
      basis[iQuad] = numpy.array([self.N1(q), self.N0(q), 
                                  self.N2(q), self.N3(q)],
                                 dtype=numpy.float64).reshape( (4,) )
      deriv = numpy.array([[self.N1p(q), self.N1q(q), self.N1r(q)],
                           [self.N0p(q), self.N0q(q), self.N0r(q)],
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

    from FIAT.reference_element import default_simplex

    cell.cellDim = 1
    shape = cell._getShape()
    self.assertEqual(default_simplex(1).get_shape(), shape.get_shape())

    cell.cellDim = 2
    shape = cell._getShape()
    self.assertEqual(default_simplex(2).get_shape(), shape.get_shape())

    cell.cellDim = 3
    shape = cell._getShape()
    self.assertEqual(default_simplex(3).get_shape(), shape.get_shape())
    return


  def test_initialize_tri3(self):
    """
    Test initialize() with tri3 cell.
    """
    cell = FIATSimplex()
    cell.inventory.dimension = 2
    cell.inventory.degree = 1
    cell._configure()
    cell.initialize(spaceDim=2)

    cellE = Tri3()
    self._checkVals(cellE, cell)
    from pylith.feassemble.CellGeometry import GeometryTri2D
    self.failUnless(isinstance(cell.geometry, GeometryTri2D))
    return


  def test_initialize_tri3_collocated(self):
    """
    Test initialize() with tri3 cell.
    """
    cell = FIATSimplex()
    cell.inventory.dimension = 2
    cell.inventory.degree = 1
    cell.inventory.collocateQuad = True
    cell._configure()
    cell.initialize(spaceDim=2)

    cellE = Tri3Collocated()
    self._checkVals(cellE, cell)
    from pylith.feassemble.CellGeometry import GeometryTri2D
    self.failUnless(isinstance(cell.geometry, GeometryTri2D))
    return


  def test_initialize_tri6(self):
    """
    Test initialize() with tri6 cell.
    """
    cell = FIATSimplex()
    cell.inventory.dimension = 2
    cell.inventory.degree = 2
    cell._configure()
    cell.initialize(spaceDim=2)

    cellE = Tri6()
    self._checkVals(cellE, cell)
    from pylith.feassemble.CellGeometry import GeometryTri2D
    self.failUnless(isinstance(cell.geometry, GeometryTri2D))
    return


  def test_initialize_tet4(self):
    """
    Test initialize() with tet4 cell.
    """
    cell = FIATSimplex()
    cell.inventory.dimension = 3
    cell.inventory.degree = 1
    cell._configure()
    cell.initialize(spaceDim=3)

    cellE = Tet4()
    self._checkVals(cellE, cell)
    from pylith.feassemble.CellGeometry import GeometryTet3D
    self.failUnless(isinstance(cell.geometry, GeometryTet3D))
    return


  def test_initialize_tet4_collocated(self):
    """
    Test initialize() with tet4 cell.
    """
    cell = FIATSimplex()
    cell.inventory.dimension = 3
    cell.inventory.degree = 1
    cell.inventory.collocateQuad = True
    cell._configure()
    cell.initialize(spaceDim=3)

    cellE = Tet4Collocated()
    self._checkVals(cellE, cell)
    from pylith.feassemble.CellGeometry import GeometryTet3D
    self.failUnless(isinstance(cell.geometry, GeometryTet3D))
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.feassemble.FIATSimplex import reference_cell
    c = reference_cell()
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
    test_scalararray(self, cellE.vertices, cell.vertices)
    test_scalararray(self, cellE.quadPts, cell.quadPts)
    test_scalararray(self, cellE.quadWts, cell.quadWts)
    test_scalararray(self, cellE.basis, cell.basis)
    test_scalararray(self, cellE.basisDeriv, cell.basisDeriv)
    return


# End of file 
