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
  return -0.25

def N0q(p):
  return -0.25

def N1(p):
  return (1+p[0])*(1-p[1])/4.0

def N1p(p):
  return +0.25

def N1q(p):
  return -0.25

def N2(p):
  return (1+p[0])*(1+p[1])/4.0

def N2p(p):
  return +0.25

def N2q(p):
  return +0.25

def N3(p):
  return (1-p[0])*(1+p[1])/4.0

def N3p(p):
  return -0.25

def N3q(p):
  return +0.25

# ----------------------------------------------------------------------
class TestFIATLagrange(unittest.TestCase):
  """
  Unit testing of FIATLagrange object.
  """
  

  def test_initialize(self):
    """
    Test initialize().
    """
    from pylith.feassemble.FIATLagrange import FIATLagrange
    cell = FIATLagrange()
    cell.degree = 1
    cell.order  = 1
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


# End of file 
