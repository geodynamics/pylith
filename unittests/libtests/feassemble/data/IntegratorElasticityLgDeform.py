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

## @file unittests/libtests/feassemble/data/IntegratorElasticityLgDeform.py

## @brief Python application for generating C++ data files for testing
## C++ elasticity integrator objects.

from IntegratorElasticity import IntegratorElasticity

import numpy
import feutils

# ----------------------------------------------------------------------

# IntegratorElasticityLgDeform class
class IntegratorElasticityLgDeform(IntegratorElasticity):
  """
  Python application for generating C++ data files for testing C++
  elasticity integrator objects.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="integratorelasticitylgdeform"):
    """
    Constructor.
    """
    IntegratorElasticity.__init__(self, name)

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
      (jacobian, jacobianInv, jacobianDet, basisDeriv) = \
                 feutils.calculateJacobian(self.quadrature, vertices)
      for iQuad in xrange(self.numQuadPts):
        wt = self.quadWts[iQuad] * jacobianDet[iQuad]
        BL0 = self._calculateBasisDerivMatLinear0(basisDeriv, iQuad)
        #cellK[:] += wt * numpy.dot(numpy.dot(BL0.transpose(), D), BL0)
        #BL1 = self._calculateBasisDerivMatLinear1(basisDeriv, iQuad)
        #cellK[:] += wt * numpy.dot(numpy.dot(BL1.transpose(), D), BL1)
        BNL = self._calculateBasisDerivMatNonlinear(basisDeriv, iQuad)
        print "BNL",BNL
        #cellK[:] += wt * numpy.dot(numpy.dot(BNL.transpose(), S), BNL)
      feutils.assembleMat(K, cellK, cell, self.spaceDim)
    return K


  def _calculateBasisDerivMatLinear0(self, basisDeriv, iQuad):
    """
    Calculate matrix of derivatives of basis functions.
    """
    if 3 == self.spaceDim:
      B = numpy.zeros( (6, self.spaceDim*self.numBasis),
                       dtype=numpy.float64)
      for iBasis in xrange(self.numBasis):
        B[0, iBasis*self.spaceDim+0] = basisDeriv[iQuad, iBasis, 0]
        B[1, iBasis*self.spaceDim+1] = basisDeriv[iQuad, iBasis, 1]
        B[2, iBasis*self.spaceDim+2] = basisDeriv[iQuad, iBasis, 2]
        B[3, iBasis*self.spaceDim+0] = basisDeriv[iQuad, iBasis, 1]
        B[3, iBasis*self.spaceDim+1] = basisDeriv[iQuad, iBasis, 0]
        B[4, iBasis*self.spaceDim+1] = basisDeriv[iQuad, iBasis, 2]
        B[4, iBasis*self.spaceDim+2] = basisDeriv[iQuad, iBasis, 1]
        B[5, iBasis*self.spaceDim+0] = basisDeriv[iQuad, iBasis, 2]
        B[5, iBasis*self.spaceDim+2] = basisDeriv[iQuad, iBasis, 0]
    elif 2 == self.spaceDim:
      B = numpy.zeros( (3, self.spaceDim*self.numBasis),
                       dtype=numpy.float64)
      for iBasis in xrange(self.numBasis):
        B[0, iBasis*self.spaceDim+0] = basisDeriv[iQuad, iBasis, 0]
        B[1, iBasis*self.spaceDim+1] = basisDeriv[iQuad, iBasis, 1]
        B[2, iBasis*self.spaceDim+0] = basisDeriv[iQuad, iBasis, 1]
        B[2, iBasis*self.spaceDim+1] = basisDeriv[iQuad, iBasis, 0]
    elif 1 == self.spaceDim:
      B = numpy.zeros( (1, self.spaceDim*self.numBasis),
                       dtype=numpy.float64)
      for iBasis in xrange(self.numBasis):
        B[0, iBasis*self.spaceDim+0] = basisDeriv[iQuad, iBasis, 0]
    else:
      raise ValueError("Unknown spatial dimension '%d'." % self.spaceDim)
    return B


  def _calculateBasisDerivMatLinear1(self, basisDeriv, iQuad):
    """
    Calculate matrix of derivatives of basis functions.
    """
    if 3 == self.spaceDim:
      B = numpy.zeros( (6, self.spaceDim*self.numBasis),
                       dtype=numpy.float64)
      for iBasis in xrange(self.numBasis):
        B[0, iBasis*self.spaceDim+0] = basisDeriv[iQuad, iBasis, 0]
        B[1, iBasis*self.spaceDim+1] = basisDeriv[iQuad, iBasis, 1]
        B[2, iBasis*self.spaceDim+2] = basisDeriv[iQuad, iBasis, 2]
        B[3, iBasis*self.spaceDim+0] = basisDeriv[iQuad, iBasis, 1]
        B[3, iBasis*self.spaceDim+1] = basisDeriv[iQuad, iBasis, 0]
        B[4, iBasis*self.spaceDim+1] = basisDeriv[iQuad, iBasis, 2]
        B[4, iBasis*self.spaceDim+2] = basisDeriv[iQuad, iBasis, 1]
        B[5, iBasis*self.spaceDim+0] = basisDeriv[iQuad, iBasis, 2]
        B[5, iBasis*self.spaceDim+2] = basisDeriv[iQuad, iBasis, 0]
    elif 2 == self.spaceDim:
      B = numpy.zeros( (3, self.spaceDim*self.numBasis),
                       dtype=numpy.float64)
      for iBasis in xrange(self.numBasis):
        B[0, iBasis*self.spaceDim+0] = basisDeriv[iQuad, iBasis, 0]
        B[1, iBasis*self.spaceDim+1] = basisDeriv[iQuad, iBasis, 1]
        B[2, iBasis*self.spaceDim+0] = basisDeriv[iQuad, iBasis, 1]
        B[2, iBasis*self.spaceDim+1] = basisDeriv[iQuad, iBasis, 0]
    elif 1 == self.spaceDim:
      B = numpy.zeros( (1, self.spaceDim*self.numBasis),
                       dtype=numpy.float64)
      for iBasis in xrange(self.numBasis):
        B[0, iBasis*self.spaceDim+0] = basisDeriv[iQuad, iBasis, 0]
    else:
      raise ValueError("Unknown spatial dimension '%d'." % self.spaceDim)
    return B


  def _calculateBasisDerivMatNonlinear(self, basisDeriv, iQuad):
    """
    Calculate matrix of derivatives of basis functions.
    """
    B = numpy.zeros( (self.spaceDim*self.spaceDim,
                      self.spaceDim*self.numBasis),
                     dtype=numpy.float64)
    for iBasis in xrange(self.numBasis):
      for iDim in xrange(self.spaceDim):
        for jDim in xrange(self.spaceDim):
          B[jDim+iDim*spaceDim, iBasis*self.spaceDim+iDim] = \
              basisDeriv[iQuad, iBasis, jDim]
    return B


# End of file 
