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

from pylith.feassemble.quadrature.Quadrature1D import Quadrature1D
from pylith.feassemble.quadrature.Quadrature1Din2D import Quadrature1Din2D
from pylith.feassemble.quadrature.Quadrature1Din3D import Quadrature1Din3D
from pylith.feassemble.quadrature.Quadrature2D import Quadrature2D
from pylith.feassemble.quadrature.Quadrature2Din3D import Quadrature2Din3D
from pylith.feassemble.quadrature.Quadrature3D import Quadrature3D


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

def verticesRef():
  return [-1.0, 1.0, 0.0]

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

    # Compute basis functions and derivatives at vertices
    basisVertE = numpy.zeros( (3, 3), dtype=numpy.float64)
    basisDerivVertE = numpy.zeros( (3, 3, 1), dtype=numpy.float64)
    iVertex = 0
    for v in verticesRef():
      basisVertE[iVertex] = numpy.array([N0(v), N1(v), N2(v)],
                                        dtype=numpy.float64).reshape( (3,) )
      deriv = numpy.array([[N0p(v)], [N1p(v)], [N2p(v)]],
                          dtype=numpy.float64)      
      basisDerivVertE[iVertex] = deriv.reshape((3, 1))
      iVertex += 1

    # Compute basis functions and derivatives at quadrature points
    basisQuadE = numpy.zeros( (2, 3), dtype=numpy.float64)
    basisDerivQuadE = numpy.zeros( (2, 3, 1), dtype=numpy.float64)
    iQuad = 0
    for q in quadPtsE:
      basisQuadE[iQuad] = numpy.array([N0(q), N1(q), N2(q)],
                                      dtype=numpy.float64).reshape( (3,) )
      deriv = numpy.array([[N0p(q)], [N1p(q)], [N2p(q)]],
                          dtype=numpy.float64)      
      basisDerivQuadE[iQuad] = deriv.reshape((3, 1))
      iQuad += 1

    quadrature = Quadrature1D()
    quadrature.cell = cell

    quadrature.initialize()

    from pylith.utils.testarray import test_double
    import pylith.feassemble.testfeassemble as testmodule

    basisVert = testmodule.basisVert(quadrature.cppHandle)
    test_double(self, basisVertE, basisVert)

    basisDerivVert = testmodule.basisDerivVert(quadrature.cppHandle)
    test_double(self, basisDerivVertE, basisDerivVert)

    basisQuad = testmodule.basisQuad(quadrature.cppHandle)
    test_double(self, basisQuadE, basisQuad)

    basisDerivQuad = testmodule.basisDerivQuad(quadrature.cppHandle)
    test_double(self, basisDerivQuadE, basisDerivQuad)

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
