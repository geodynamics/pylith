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

## @file unittests/pytests/feassemble/TestQuadrature.py

## @brief Unit testing of Python Quadrature object.

import unittest
import numpy
from pylith.feassemble.Quadrature import *

# ----------------------------------------------------------------------
def N0(p):
  return -0.5*p*(1.0-p)

def N0p(p):
  return -0.5*(1.0-p) + 0.5*p

def N1(p):
  return 0.5*p*(1.0+p)

def N1p(p):
  return +0.5*(1.0+p) + 0.5*p

def N2(p):
  return (1.0-p**2)

def N2p(p):
  return -2.0*p


# ----------------------------------------------------------------------
class TestQuadrature(unittest.TestCase):
  """
  Unit testing of Python Quadrature object.
  """

  def test_minJacobian(self):
    """
    Test minJacobian attribute.
    """
    minJacobian = 4.0e-02;
    q = Quadrature1D()
    q.minJacobian = minJacobian
    self.assertEqual(minJacobian, q.minJacobian)
    return
    

  def test_initialize(self):
    """
    Test initialize().
    """
    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATSimplex()
    cell.shape = "line"
    cell.degree = 2
    cell.order = 2
    quadPtsE = numpy.array( [(-1.0/3**0.5,),
                             (+1.0/3**0.5,)],
                            dtype=numpy.float64 )
    quadWtsE = numpy.array( [1.0, 1.0], dtype=numpy.float64 )
    basisE = numpy.zeros( (2, 3), dtype=numpy.Float64)
    basisDerivE = numpy.zeros( (2, 3, 1), dtype=numpy.Float64)
    iQuad = 0
    for q in quadPtsE:
      basisE[iQuad] = numpy.array([N0(q), N1(q), N2(q)],
                                 dtype=numpy.float64).reshape( (3,) )
      deriv = numpy.array([[N0p(q)], [N1p(q)], [N2p(q)]],
                          dtype=numpy.Float64)      
      basisDerivE[iQuad] = deriv.reshape((3, 1))
      iQuad += 1

    quadrature = Quadrature1D()
    quadrature.cell = cell

    quadrature.initialize()

    from pylith.utils.testarray import test_double
    import pylith.feassemble.testfeassemble as testmodule

    basis = testmodule.basis(quadrature.cppHandle)
    test_double(self, basisE, basis)

    basisDeriv = testmodule.basisDeriv(quadrature.cppHandle)
    test_double(self, basisDerivE, basisDeriv)

    quadWts = testmodule.quadWts(quadrature.cppHandle)
    test_double(self, quadWtsE, quadWts)
    
    quadPts = testmodule.quadPtsRef(quadrature.cppHandle)
    test_double(self, quadPtsE, quadPts)
    
    return


  def test_constructors(self):
    """
    Test constructors for quadrature objects.
    """
    q = Quadrature1D()
    self.assertEqual(1, q.spaceDim)
    self.assertEqual(1, q.cellDim)
    self.failIf(None == q.cppHandle)
    
    q = Quadrature1Din2D()
    self.assertEqual(2, q.spaceDim)
    self.assertEqual(1, q.cellDim)
    self.failIf(None == q.cppHandle)
    
    q = Quadrature1Din3D()
    self.assertEqual(3, q.spaceDim)
    self.assertEqual(1, q.cellDim)
    self.failIf(None == q.cppHandle)
    
    q = Quadrature2D()
    self.assertEqual(2, q.spaceDim)
    self.assertEqual(2, q.cellDim)
    self.failIf(None == q.cppHandle)
    
    q = Quadrature2Din3D()
    self.assertEqual(3, q.spaceDim)
    self.assertEqual(2, q.cellDim)
    self.failIf(None == q.cppHandle)
    
    q = Quadrature3D()
    self.assertEqual(3, q.spaceDim)
    self.assertEqual(3, q.cellDim)
    self.failIf(None == q.cppHandle)
    
    return


# End of file 
