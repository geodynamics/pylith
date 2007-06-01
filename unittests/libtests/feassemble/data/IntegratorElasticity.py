#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

## @file unittests/libtests/feassemble/data/IntegratorElasticity.py

## @brief Python application for generating C++ data files for testing
## C++ elasticity integrator objects.

from IntegratorApp import IntegratorApp

import numpy
import feutils

# ----------------------------------------------------------------------

# IntegratorElasticity class
class IntegratorElasticity(IntegratorApp):
  """
  Python application for generating C++ data files for testing C++
  elasticity integrator objects.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="integratorelasticity"):
    """
    Constructor.
    """
    IntegratorApp.__init__(self, name)

    self.density = None
    self.lameMu = None
    self.lameLambda = None

    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _calculateStiffnessMat(self):
    """
    Calculate stiffness matrix.

    """
    import feutils

    K = numpy.zeros( (self.spaceDim*self.numVertices,
                      self.spaceDim*self.numVertices),
                     dtype=numpy.float64)

    # Matrix of elasticity values
    D = self._calculateElasticityMat()
    
    for cell in self.cells:
      cellK = numpy.zeros( (self.spaceDim*self.numBasis,
                            self.spaceDim*self.numBasis),
                           dtype=numpy.float64)
      vertices = self.vertices[cell, :]
      (jacobian, jacobianInv, jacobianDet) = \
                 feutils.calculateJacobian(self.quadrature, vertices)
      for iQuad in xrange(self.numQuadPts):
        wt = self.quadWts[iQuad] * jacobianDet[iQuad]
        B = self._calculateBasisDerivMat(iQuad)
        cellK[:] += wt * numpy.dot(numpy.dot(B.transpose(), D), B)
      feutils.assembleMat(K, cellK, cell, self.spaceDim)
    return K


  def _calculateMassMat(self):
    """
    Calculate mass matrix.
    """

    M = numpy.zeros( (self.spaceDim*self.numVertices,
                      self.spaceDim*self.numVertices),
                     dtype=numpy.float64)

    for cell in self.cells:
      cellM = numpy.zeros( (self.spaceDim*self.numBasis,
                            self.spaceDim*self.numBasis),
                           dtype=numpy.float64)
      vertices = self.vertices[cell, :]
      (jacobian, jacobianInv, jacobianDet) = \
                 feutils.calculateJacobian(self.quadrature, vertices)
      for iQuad in xrange(self.numQuadPts):
        wt = self.quadWts[iQuad] * jacobianDet[iQuad]
        N = self._calculateBasisMat(iQuad)
        cellM[:] += self.density * wt * numpy.dot(N.transpose(), N)
      feutils.assembleMat(M, cellM, cell, self.spaceDim)
    return M


  def _calculateElasticityMat(self):
    """
    Calculate elasticity matrix.
    """
    if 1 == self.cellDim:
      lambda2mu = self.lameLambda + 2*self.lameMu
      C1111 = lambda2mu
      D = numpy.array([ [C1111] ],
                      dtype=numpy.float64)
    elif 2 == self.cellDim:
      lambda2mu = self.lameLambda + 2*self.lameMu
      C1111 = lambda2mu
      C1122 = self.lameLambda
      C1112 = 0.0
      C2222 = lambda2mu
      C2212 = 0.0
      C1212 = 2.0*self.lameMu
      D = numpy.array([ [C1111, C1122, C1112],
                        [C1122, C2222, C2212],
                        [C1112, C2212, C1212] ],
                      dtype=numpy.float64)
    elif 3 == self.cellDim:
      lambda2mu = self.lameLambda + 2.0*self.lameMu
      C1111 = lambda2mu
      C1122 = self.lameLambda
      C1133 = self.lameLambda
      C1112 = 0.0
      C1123 = 0.0
      C1113 = 0.0
      C2222 = lambda2mu
      C2233 = self.lameLambda
      C2212 = 0.0
      C2223 = 0.0
      C2213 = 0.0
      C3333 = lambda2mu
      C3312 = 0.0
      C3323 = 0.0
      C3313 = 0.0
      C1212 = 2.0*self.lameMu
      C1223 = 0.0
      C1213 = 0.0
      C2323 = 2.0*self.lameMu
      C2313 = 0.0
      C1313 = 2.0*self.lameMu
      D = numpy.array([ [C1111, C1122, C1133, C1112, C1123, C1113],
                        [C1122, C2222, C2233, C2212, C2223, C2213],
                        [C1133, C2233, C3333, C3312, C3323, C3313],
                        [C1112, C2212, C3312, C1212, C1223, C1213],
                        [C1123, C2223, C3323, C1223, C2323, C2313],
                        [C1113, C2213, C3313, C1213, C2313, C1313] ],
                      dtype=numpy.float64)
    return D


  def _calculateBasisMat(self, iQuad):
    """
    Calculate matrix of basis functions.
    """
    N = numpy.zeros( (self.spaceDim, self.spaceDim*self.numBasis),
                     dtype=numpy.float64)
    for iBasis in xrange(self.numBasis):
      for iDim in xrange(self.spaceDim):
        N[iDim, iBasis*self.spaceDim+iDim] = self.basis[iQuad, iBasis]
    return N


  def _calculateBasisDerivMat(self, iQuad):
    """
    Calculate matrix of derivatives of basis functions.
    """
    if 3 == self.spaceDim:
      B = numpy.zeros( (6, self.spaceDim*self.numBasis),
                       dtype=numpy.float64)
      for iBasis in xrange(self.numBasis):
        B[0, iBasis*self.spaceDim+0] = self.basisDeriv[iQuad, iBasis, 0]
        B[1, iBasis*self.spaceDim+1] = self.basisDeriv[iQuad, iBasis, 1]
        B[2, iBasis*self.spaceDim+2] = self.basisDeriv[iQuad, iBasis, 2]
        B[3, iBasis*self.spaceDim+0] = 0.5*self.basisDeriv[iQuad, iBasis, 1]
        B[3, iBasis*self.spaceDim+1] = 0.5*self.basisDeriv[iQuad, iBasis, 0]
        B[4, iBasis*self.spaceDim+1] = 0.5*self.basisDeriv[iQuad, iBasis, 2]
        B[4, iBasis*self.spaceDim+2] = 0.5*self.basisDeriv[iQuad, iBasis, 1]
        B[5, iBasis*self.spaceDim+0] = 0.5*self.basisDeriv[iQuad, iBasis, 2]
        B[5, iBasis*self.spaceDim+2] = 0.5*self.basisDeriv[iQuad, iBasis, 0]
    elif 2 == self.spaceDim:
      B = numpy.zeros( (3, self.spaceDim*self.numBasis),
                       dtype=numpy.float64)
      for iBasis in xrange(self.numBasis):
        B[0, iBasis*self.spaceDim+0] = self.basisDeriv[iQuad, iBasis, 0]
        B[1, iBasis*self.spaceDim+1] = self.basisDeriv[iQuad, iBasis, 1]
        B[2, iBasis*self.spaceDim+0] = 0.5*self.basisDeriv[iQuad, iBasis, 1]
        B[2, iBasis*self.spaceDim+1] = 0.5*self.basisDeriv[iQuad, iBasis, 0]
    elif 1 == self.spaceDim:
      B = numpy.zeros( (1, self.spaceDim*self.numBasis),
                       dtype=numpy.float64)
      for iBasis in xrange(self.numBasis):
        B[0, iBasis*self.spaceDim+0] = self.basisDeriv[iQuad, iBasis, 0]
    else:
      raise ValueError("Unknown spatial dimension '%d'." % self.spaceDim)
    return B


# End of file 
