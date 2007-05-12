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
from pylith.utils.testarray import test_double

# ----------------------------------------------------------------------
def N0(p):
  return -0.5*p*(1.0-p)

def N0p(p):
  return 1.0*p - 0.5

def N1(p):
  return 0.5*p*(1.0+p)

def N1p(p):
  return 1.0*p + 0.5

def N2(p):
  return (1.0-p*p)

def N2p(p):
  return -2.0*p

def nodesRef():
  return [-1.0, 1.0, 0.0]


# ----------------------------------------------------------------------
class TestFIATSimplex(unittest.TestCase):
  """
  Unit testing of FIATSimplex object.
  """
  

  def test_shape(self):
    """
    Test _getShape().
    """
    from pylith.feassemble.FIATSimplex import FIATSimplex
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


  def test_initialize(self):
    """
    Test initialize().
    """
    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATSimplex()
    cell.shape  = "line"
    cell.degree = 2
    cell.order  = 2
    quadPtsE = numpy.array( [(-1.0/3**0.5,),
                             (+1.0/3**0.5,)],
                            dtype=numpy.float64 )
    quadWtsE = numpy.array( [1.0, 1.0], dtype=numpy.float64 )

    # Compute basis fns and derivatives at quadrature points
    basis = numpy.zeros( (2, 3), dtype=numpy.float64)
    basisDeriv = numpy.zeros( (2, 3, 1), dtype=numpy.float64)
    iQuad = 0
    for q in quadPtsE:
      basis[iQuad] = numpy.array([N0(q), N1(q), N2(q)],
                                 dtype=numpy.float64).reshape( (3,) )
      deriv = numpy.array([[N0p(q)], [N1p(q)], [N2p(q)]],
                          dtype=numpy.float64)      
      basisDeriv[iQuad] = deriv.reshape((3, 1))
      iQuad += 1

    cell.initialize()

    # Check basic attributes
    self.assertEqual(1, cell.cellDim)
    self.assertEqual(3, cell.numCorners)
    self.assertEqual(2, cell.numQuadPts)

    # Check arrays
    test_double(self, quadPtsE, cell.quadPts)
    test_double(self, quadWtsE, cell.quadWts)
    test_double(self, basis, cell.basis)
    test_double(self, basisDeriv, cell.basisDeriv)

    return


# End of file 
