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
from pylith.utils.testarray import test_double

# ----------------------------------------------------------------------
def N0(p):
  return (1-p[0])*(1-p[1])/4.0

def N0p(p):
  return -(1-p[1])/4.0

def N0q(p):
  return -(1-p[0])/4.0

def N1(p):
  return (1+p[0])*(1-p[1])/4.0

def N1p(p):
  return (1-p[1])/4.0

def N1q(p):
  return -(1+p[0])/4.0

def N2(p):
  return (1-p[0])*(1+p[1])/4.0

def N2p(p):
  return -(1+p[1])/4.0

def N2q(p):
  return (1-p[0])/4.0

def N3(p):
  return (1+p[0])*(1+p[1])/4.0

def N3p(p):
  return (1+p[1])/4.0

def N3q(p):
  return (1+p[0])/4.0

def Q0(p):
  return (1-p[0])*(1-p[1])*(1-p[2])/8.0

def Q0p(p):
  return -(1-p[1])*(1-p[2])/8.0

def Q0q(p):
  return -(1-p[0])*(1-p[2])/8.0

def Q0r(p):
  return -(1-p[0])*(1-p[1])/8.0

def Q1(p):
  return (1+p[0])*(1-p[1])*(1-p[2])/8.0

def Q1p(p):
  return (1-p[1])*(1-p[2])/8.0

def Q1q(p):
  return -(1+p[0])*(1-p[2])/8.0

def Q1r(p):
  return -(1+p[0])*(1-p[1])/8.0

def Q2(p):
  return (1-p[0])*(1+p[1])*(1-p[2])/8.0

def Q2p(p):
  return -(1+p[1])*(1-p[2])/8.0

def Q2q(p):
  return (1-p[0])*(1-p[2])/8.0

def Q2r(p):
  return -(1-p[0])*(1+p[1])/8.0

def Q3(p):
  return (1+p[0])*(1+p[1])*(1-p[2])/8.0

def Q3p(p):
  return (1+p[1])*(1-p[2])/8.0

def Q3q(p):
  return (1+p[0])*(1-p[2])/8.0

def Q3r(p):
  return -(1+p[0])*(1+p[1])/8.0

def Q4(p):
  return (1-p[0])*(1-p[1])*(1+p[2])/8.0

def Q4p(p):
  return -(1-p[1])*(1+p[2])/8.0

def Q4q(p):
  return -(1-p[0])*(1+p[2])/8.0

def Q4r(p):
  return (1-p[0])*(1-p[1])/8.0

def Q5(p):
  return (1+p[0])*(1-p[1])*(1+p[2])/8.0

def Q5p(p):
  return (1-p[1])*(1+p[2])/8.0

def Q5q(p):
  return -(1+p[0])*(1+p[2])/8.0

def Q5r(p):
  return (1+p[0])*(1-p[1])/8.0

def Q6(p):
  return (1-p[0])*(1+p[1])*(1+p[2])/8.0

def Q6p(p):
  return -(1+p[1])*(1+p[2])/8.0

def Q6q(p):
  return (1-p[0])*(1+p[2])/8.0

def Q6r(p):
  return (1-p[0])*(1+p[1])/8.0

def Q7(p):
  return (1+p[0])*(1+p[1])*(1+p[2])/8.0

def Q7p(p):
  return (1+p[1])*(1+p[2])/8.0

def Q7q(p):
  return (1+p[0])*(1+p[2])/8.0

def Q7r(p):
  return (1+p[0])*(1+p[1])/8.0

# ----------------------------------------------------------------------
class TestFIATLagrange(unittest.TestCase):
  """
  Unit testing of FIATLagrange object.
  """


  def test_initialize_quad(self):
    """
    Test initialize().
    """
    from pylith.feassemble.FIATLagrange import FIATLagrange
    cell = FIATLagrange()
    cell.cellDim = 2
    cell.degree = 1
    cell.order  = 2
    quadPtsE = numpy.array( [(-1.0/3**0.5, -1.0/3**0.5),
                             (+1.0/3**0.5, -1.0/3**0.5),
                             (-1.0/3**0.5, +1.0/3**0.5),
                             (+1.0/3**0.5, +1.0/3**0.5)],
                            dtype=numpy.float64 )
    quadWtsE = numpy.array( [1.0, 1.0, 1.0, 1.0], dtype=numpy.float64 )

    # Compute basis fns and derivatives at quadrature points
    basis = numpy.zeros( (4, 4), dtype=numpy.float64)
    basisDeriv = numpy.zeros( (4, 4, 2), dtype=numpy.float64)
    iQuad = 0
    for q in quadPtsE:
      basis[iQuad] = numpy.array([N0(q), N1(q), N2(q), N3(q)],
                                 dtype=numpy.float64).reshape( (4,) )
      deriv = numpy.array([[N0p(q), N0q(q)],
                           [N1p(q), N1q(q)],
                           [N2p(q), N2q(q)],
                           [N3p(q), N3q(q)]],
                          dtype=numpy.float64)
      basisDeriv[iQuad] = deriv.reshape((4, 2))
      iQuad += 1

    cell.initialize()

    # Check basic attributes
    self.assertEqual(2, cell.cellDim)
    self.assertEqual(4, cell.numCorners)
    self.assertEqual(4, cell.numQuadPts)

    # Check arrays
    test_double(self, quadPtsE, cell.quadPts)
    test_double(self, quadWtsE, cell.quadWts)
    test_double(self, basis, cell.basis)
    test_double(self, basisDeriv, cell.basisDeriv)

    return


  def test_initialize_hex(self):
    """
    Test initialize().
    """
    from pylith.feassemble.FIATLagrange import FIATLagrange
    cell = FIATLagrange()
    cell.cellDim = 3
    cell.degree = 1
    cell.order  = 2
    quadPtsE = numpy.array( [(-1.0/3**0.5, -1.0/3**0.5, -1.0/3**0.5),
                             (+1.0/3**0.5, -1.0/3**0.5, -1.0/3**0.5),
                             (-1.0/3**0.5, +1.0/3**0.5, -1.0/3**0.5),
                             (+1.0/3**0.5, +1.0/3**0.5, -1.0/3**0.5),
                             (-1.0/3**0.5, -1.0/3**0.5, +1.0/3**0.5),
                             (+1.0/3**0.5, -1.0/3**0.5, +1.0/3**0.5),
                             (-1.0/3**0.5, +1.0/3**0.5, +1.0/3**0.5),
                             (+1.0/3**0.5, +1.0/3**0.5, +1.0/3**0.5)],
                            dtype=numpy.float64 )
    quadWtsE = numpy.array( [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], dtype=numpy.float64 )

    # Compute basis fns and derivatives at quadrature points
    iQuad = 0
    basis = numpy.zeros( (8, 8), dtype=numpy.float64)
    basisDeriv = numpy.zeros( (8, 8, 3), dtype=numpy.float64)
    for q in quadPtsE:
      basis[iQuad] = numpy.array([Q0(q), Q1(q), Q2(q), Q3(q), Q4(q), Q5(q), Q6(q), Q7(q)],
                                 dtype=numpy.float64).reshape( (8,) )
      deriv = numpy.array([[Q0p(q), Q0q(q), Q0r(q)],
                           [Q1p(q), Q1q(q), Q1r(q)],
                           [Q2p(q), Q2q(q), Q2r(q)],
                           [Q3p(q), Q3q(q), Q3r(q)],
                           [Q4p(q), Q4q(q), Q4r(q)],
                           [Q5p(q), Q5q(q), Q5r(q)],
                           [Q6p(q), Q6q(q), Q6r(q)],
                           [Q7p(q), Q7q(q), Q7r(q)]],
                          dtype=numpy.float64)
      basisDeriv[iQuad] = deriv.reshape((8, 3))
      iQuad += 1

    cell.initialize()

    # Check basic attributes
    self.assertEqual(3, cell.cellDim)
    self.assertEqual(8, cell.numCorners)
    self.assertEqual(8, cell.numQuadPts)

    # Check arrays
    test_double(self, quadPtsE, cell.quadPts)
    test_double(self, quadWtsE, cell.quadWts)
    test_double(self, basis, cell.basis)
    test_double(self, basisDeriv, cell.basisDeriv)

    return


# End of file 
