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
      fieldTpdt = self.fieldT + self.fieldTIncr
      for iQuad in xrange(self.numQuadPts):
        wt = self.quadWts[iQuad] * jacobianDet[iQuad]
        BL0 = self._calculateBasisDerivMatLinear0(basisDeriv, iQuad)
        cellK[:] += wt * numpy.dot(numpy.dot(BL0.transpose(), D), BL0)
        BL1 = self._calculateBasisDerivMatLinear1(basisDeriv, iQuad, fieldTpdt)
        cellK[:] += wt * numpy.dot(numpy.dot(BL1.transpose(), D), BL1)
        BNL = self._calculateBasisDerivMatNonlinear(basisDeriv, iQuad)
        strain = numpy.dot(BL0+BL1, fieldTpdt)
        S = self._calculateStress(strain, D)
        cellK[:] += wt * numpy.dot(numpy.dot(BNL.transpose(), S), BNL)
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


  def _calculateBasisDerivMatLinear1(self, basisDeriv, iQuad, disp):
    """
    Calculate matrix of derivatives of basis functions.
    """
    if 3 == self.spaceDim:
      B = numpy.zeros( (6, self.spaceDim*self.numBasis),
                       dtype=numpy.float64)
      l11 = 0.0
      l12 = 0.0
      l13 = 0.0
      l21 = 0.0
      l22 = 0.0
      l23 = 0.0
      l31 = 0.0
      l32 = 0.0
      l33 = 0.0
      for kBasis in xrange(self.numBasis):
        l11 += basisDeriv[iQuad, kBasis, 0]*disp[kBasis*self.spaceDim  ]
        l12 += basisDeriv[iQuad, kBasis, 1]*disp[kBasis*self.spaceDim  ]
        l13 += basisDeriv[iQuad, kBasis, 2]*disp[kBasis*self.spaceDim  ]
        l21 += basisDeriv[iQuad, kBasis, 0]*disp[kBasis*self.spaceDim+1]
        l22 += basisDeriv[iQuad, kBasis, 1]*disp[kBasis*self.spaceDim+1]
        l23 += basisDeriv[iQuad, kBasis, 2]*disp[kBasis*self.spaceDim+1]
        l31 += basisDeriv[iQuad, kBasis, 0]*disp[kBasis*self.spaceDim+2]
        l32 += basisDeriv[iQuad, kBasis, 1]*disp[kBasis*self.spaceDim+2]
        l33 += basisDeriv[iQuad, kBasis, 2]*disp[kBasis*self.spaceDim+2]
      for iBasis in xrange(self.numBasis):
        B[0, iBasis*self.spaceDim+0] = basisDeriv[iQuad, iBasis, 0]*l11
        B[0, iBasis*self.spaceDim+1] = basisDeriv[iQuad, iBasis, 0]*l21
        B[0, iBasis*self.spaceDim+2] = basisDeriv[iQuad, iBasis, 0]*l31
        B[1, iBasis*self.spaceDim+0] = basisDeriv[iQuad, iBasis, 1]*l12
        B[1, iBasis*self.spaceDim+1] = basisDeriv[iQuad, iBasis, 1]*l22
        B[1, iBasis*self.spaceDim+2] = basisDeriv[iQuad, iBasis, 1]*l32
        B[2, iBasis*self.spaceDim+0] = basisDeriv[iQuad, iBasis, 2]*l13
        B[2, iBasis*self.spaceDim+1] = basisDeriv[iQuad, iBasis, 2]*l23
        B[2, iBasis*self.spaceDim+2] = basisDeriv[iQuad, iBasis, 2]*l33

        B[3, iBasis*self.spaceDim+0] = \
            basisDeriv[iQuad, iBasis, 1]*l11 + basisDeriv[iQuad, iBasis, 0]*l12
        B[3, iBasis*self.spaceDim+1] = \
            basisDeriv[iQuad, iBasis, 0]*l22 + basisDeriv[iQuad, iBasis, 1]*l21
        B[3, iBasis*self.spaceDim+2] = \
            basisDeriv[iQuad, iBasis, 1]*l31 + basisDeriv[iQuad, iBasis, 0]*l32
        
        B[4, iBasis*self.spaceDim+0] = \
            basisDeriv[iQuad, iBasis, 2]*l12 + basisDeriv[iQuad, iBasis, 1]*l13
        B[4, iBasis*self.spaceDim+1] = \
            basisDeriv[iQuad, iBasis, 2]*l22 + basisDeriv[iQuad, iBasis, 1]*l23
        B[4, iBasis*self.spaceDim+2] = \
            basisDeriv[iQuad, iBasis, 1]*l33 + basisDeriv[iQuad, iBasis, 2]*l32
        
        B[5, iBasis*self.spaceDim+0] = \
            basisDeriv[iQuad, iBasis, 2]*l11 + basisDeriv[iQuad, iBasis, 0]*l13
        B[5, iBasis*self.spaceDim+1] = \
            basisDeriv[iQuad, iBasis, 2]*l21 + basisDeriv[iQuad, iBasis, 0]*l23
        B[5, iBasis*self.spaceDim+2] = \
            basisDeriv[iQuad, iBasis, 0]*l33 + basisDeriv[iQuad, iBasis, 2]*l31

    elif 2 == self.spaceDim:
      B = numpy.zeros( (3, self.spaceDim*self.numBasis),
                       dtype=numpy.float64)
      l11 = 0.0
      l12 = 0.0
      l21 = 0.0
      l22 = 0.0
      for kBasis in xrange(self.numBasis):
        l11 += basisDeriv[iQuad, kBasis, 0]*disp[kBasis*self.spaceDim  ]
        l12 += basisDeriv[iQuad, kBasis, 1]*disp[kBasis*self.spaceDim  ]
        l21 += basisDeriv[iQuad, kBasis, 0]*disp[kBasis*self.spaceDim+1]
        l22 += basisDeriv[iQuad, kBasis, 1]*disp[kBasis*self.spaceDim+1]
      for iBasis in xrange(self.numBasis):
        B[0, iBasis*self.spaceDim+0] = basisDeriv[iQuad, iBasis, 0]*l11
        B[0, iBasis*self.spaceDim+1] = basisDeriv[iQuad, iBasis, 0]*l21
        B[1, iBasis*self.spaceDim+0] = basisDeriv[iQuad, iBasis, 1]*l12
        B[1, iBasis*self.spaceDim+1] = basisDeriv[iQuad, iBasis, 1]*l22
        B[2, iBasis*self.spaceDim+0] = \
            basisDeriv[iQuad, iBasis, 1]*l11 + basisDeriv[iQuad, iBasis, 0]*l12
        B[2, iBasis*self.spaceDim+1] = \
            basisDeriv[iQuad, iBasis, 0]*l22 + basisDeriv[iQuad, iBasis, 1]*l21

    elif 1 == self.spaceDim:
      B = numpy.zeros( (1, self.spaceDim*self.numBasis),
                       dtype=numpy.float64)
      l11 = 0.0
      for kBasis in xrange(self.numBasis):
        l11 += basisDeriv[iQuad, kBasis, 0]*disp[kBasis]
      for iBasis in xrange(self.numBasis):
        B[0, iBasis*self.spaceDim+0] = basisDeriv[iQuad, iBasis, 0]*l11
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
          B[jDim+iDim*self.spaceDim, iBasis*self.spaceDim+iDim] = \
              basisDeriv[iQuad, iBasis, jDim]
    return B


  def _calculateStress(self, strain, D):
    """
    Calculte 2nd Priola-Kirchoff stress matrix.
    """
    S = numpy.zeros( (self.spaceDim*self.spaceDim,
                      self.spaceDim*self.spaceDim), dtype=numpy.float64)
    Svec = numpy.dot(D, strain)
    if 3 == self.spaceDim:
      Smat = numpy.array([[Svec[0], Svec[3], Svec[5]],
                          [Svec[3], Svec[1], Svec[4]],
                          [Svec[5], Svec[4], Svec[2]]], dtype=numpy.float64)
      S[0:3,0:3] = Smat[:]
      S[3:6,3:6] = Smat[:]
      S[6:9,6:9] = Smat[:]
    elif 2 == self.spaceDim:
      Smat = numpy.array([[Svec[0], Svec[2]],
                          [Svec[2], Svec[1]]], dtype=numpy.float64)
      S[0:2,0:2] = Smat[:]
      S[2:4,2:4] = Smat[:]
    elif 1 == self.spaceDim:
      Smat = numpy.array([[Svec[0]]], dtype=numpy.float64)
      S[0:1,0:1] = Smat[:]
    return S


# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = IntegratorElasticityLgDeform()
  app.run()


# End of file 
