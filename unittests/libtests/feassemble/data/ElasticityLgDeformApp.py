#!/usr/bin/env python
#
# ----------------------------------------------------------------------
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
# ----------------------------------------------------------------------
#

## @file unittests/libtests/feassemble/data/ElasticityLgDeformApp.py

## @brief Python application for generating C++ data files for testing
## C++ elasticity integrator objects.

from ElasticityApp import ElasticityApp

import numpy
import feutils

# ----------------------------------------------------------------------

# ElasticityLgDeformApp class
class ElasticityLgDeformApp(ElasticityApp):
  """
  Python application for generating C++ data files for testing C++
  elasticity integrator objects.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="elasticitylgdeformapp"):
    """
    Constructor.
    """
    ElasticityApp.__init__(self, name)

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
        BL1 = self._calculateBasisDerivMatLinear1(basisDeriv, iQuad, fieldTpdt)
        BL = BL0 + BL1
        cellK[:] += wt * numpy.dot(numpy.dot(BL.transpose(), D), BL)
        BNL = self._calculateBasisDerivMatNonlinear(basisDeriv, iQuad)
        strain = self._calculateStrain(basisDeriv, iQuad, fieldTpdt)
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


  def _calculateStrain(self, basisDeriv, iQuad, disp):
    """
    Calculte Green-Lagrange strain. Shear strains are twice the
    Green-Lagrance values for compatibility with computing the strains
    using the B matrix in the infinitesimal strain case.
    """
    if 3 == self.spaceDim:
      strain = numpy.zeros( (1,6), dtype=numpy.float64)

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
      strain[0, 0] = 0.5*(l11*l11 + l21*l21 + l31*l31)
      strain[0, 1] = 0.5*(l12*l12 + l22*l22 + l32*l32)
      strain[0, 2] = 0.5*(l13*l13 + l23*l23 + l33*l33)
      strain[0, 3] = (l11*l12 + l21*l22 + l31*l32) # Use 2*e12 (D has mu)
      strain[0, 4] = (l12*l13 + l22*l23 + l32*l33)
      strain[0, 5] = (l11*l13 + l21*l23 + l31*l33)
      for iBasis in xrange(self.numBasis):
        strain[0, 0] += \
            basisDeriv[iQuad, iBasis, 0]*disp[iBasis*self.spaceDim  ]
        strain[0, 1] += \
            basisDeriv[iQuad, iBasis, 1]*disp[iBasis*self.spaceDim+1]
        strain[0, 2] += \
            basisDeriv[iQuad, iBasis, 2]*disp[iBasis*self.spaceDim+2]
        strain[0, 3] += \
            (basisDeriv[iQuad, iBasis, 0]*disp[iBasis*self.spaceDim+1] +
             basisDeriv[iQuad, iBasis, 1]*disp[iBasis*self.spaceDim  ])
        strain[0, 4] += \
            (basisDeriv[iQuad, iBasis, 1]*disp[iBasis*self.spaceDim+2] +
             basisDeriv[iQuad, iBasis, 2]*disp[iBasis*self.spaceDim+1])
        strain[0, 5] += \
            (basisDeriv[iQuad, iBasis, 0]*disp[iBasis*self.spaceDim+2] +
             basisDeriv[iQuad, iBasis, 2]*disp[iBasis*self.spaceDim  ])

    elif 2 == self.spaceDim:
      strain = numpy.zeros( (1,3), dtype=numpy.float64)
      l11 = 0.0
      l12 = 0.0
      l21 = 0.0
      l22 = 0.0
      for kBasis in xrange(self.numBasis):
        l11 += basisDeriv[iQuad, kBasis, 0]*disp[kBasis*self.spaceDim  ]
        l12 += basisDeriv[iQuad, kBasis, 1]*disp[kBasis*self.spaceDim  ]
        l21 += basisDeriv[iQuad, kBasis, 0]*disp[kBasis*self.spaceDim+1]
        l22 += basisDeriv[iQuad, kBasis, 1]*disp[kBasis*self.spaceDim+1]
      strain[0, 0] = 0.5*(l11*l11 + l21*l21)
      strain[0, 1] = 0.5*(l12*l12 + l22*l22)
      strain[0, 2] = (l11*l12 + l21*l22) # Use 2*e12 (D has mu, not 2*mu)
      for iBasis in xrange(self.numBasis):
        strain[0, 0] += \
            basisDeriv[iQuad, iBasis, 0]*disp[iBasis*self.spaceDim  ]
        strain[0, 1] += \
            basisDeriv[iQuad, iBasis, 1]*disp[iBasis*self.spaceDim+1]
        strain[0, 2] += \
            (basisDeriv[iQuad, iBasis, 0]*disp[iBasis*self.spaceDim+1] +
             basisDeriv[iQuad, iBasis, 1]*disp[iBasis*self.spaceDim  ])

    elif 1 == self.spaceDim:
      strain = numpy.zeros( (1,1), dtype=numpy.float64)
      l11 = 0.0
      for kBasis in xrange(self.numBasis):
        l11 += basisDeriv[iQuad, kBasis, 0]*disp[kBasis]
      strain[0, 0] = 0.5*l11*l11
      for iBasis in xrange(self.numBasis):
        strain[0, 0] += basisDeriv[iQuad, iBasis, 0]*disp[iBasis]
    else:
      raise ValueError("Unknown spatial dimension '%d'." % self.spaceDim)
      
    return strain


  def _calculateStress(self, strain, D):
    """
    Calculte 2nd Priola-Kirchoff stress matrix.
    """
    S = numpy.zeros( (self.spaceDim*self.spaceDim,
                      self.spaceDim*self.spaceDim), dtype=numpy.float64)
    Svec = numpy.dot(D, strain.transpose())
    if 3 == self.spaceDim:
      Smat = numpy.array([[Svec[0,0], Svec[3,0], Svec[5,0]],
                          [Svec[3,0], Svec[1,0], Svec[4,0]],
                          [Svec[5,0], Svec[4,0], Svec[2,0]]], 
                         dtype=numpy.float64)
      S[0:3,0:3] = Smat[:]
      S[3:6,3:6] = Smat[:]
      S[6:9,6:9] = Smat[:]
    elif 2 == self.spaceDim:
      Smat = numpy.array([[Svec[0,0], Svec[2,0]],
                          [Svec[2,0], Svec[1,0]]], dtype=numpy.float64)
      S[0:2,0:2] = Smat[:]
      S[2:4,2:4] = Smat[:]
    elif 1 == self.spaceDim:
      Smat = numpy.array([[Svec[0]]], dtype=numpy.float64)
      S[0:1,0:1] = Smat[:]
    else:
      raise ValueError("Unknown spatial dimension '%d'." % self.spaceDim)
    return S


# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = ElasticityLgDeformApp()
  app.run()


# End of file 
