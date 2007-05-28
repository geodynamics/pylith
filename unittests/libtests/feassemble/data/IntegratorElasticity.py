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

    self.K = self._calculateStiffnessMat()
    self.M = self._calculateMassMat()
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _calculateStiffnessMat(self):
    """
    Calculate stiffness matrix.

    """
    import feutils

    self.K = numpy.zeros ( (self.spaceDim*self.numVertices,
                            self.spaceDim*self.numVertices),
                           dtype=numpy.float64)

    # Matrix of elasticity values
    D = self._calculateElasticityMat()
    
    for cell in self.cells:
      vertices = self.vertices[cell, :]

      # Matrix of shape function derivatives
      B = self._calculateBasisDerivMat(vertices)

      cellK = numpy.dot(numpy.dot(B.transpose(), D), B)

      feutils.assembleMat(self.K, cellK, cell, self.spaceDim)
    return


  def _calculateMassMat(self):
    """
    Calculate mass matrix.
    """

    self.M = numpy.zeros ( (self.spaceDim*self.numVertices,
                            self.spaceDim*self.numVertices),
                           dtype=numpy.float64)


    for cell in self.cells:
      vertices = self.vertices[cell, :]
      (jacobian, jacobianInv, jacobianDet) = \
                 feutils.calculateJacobianQuad(self.quadrature, vertices)
      for iQuad in xrange(self.numQuadPts):
        wt = self.quadWts[iQuad] * jacobianDet[iQuad]
        N = self._calculateBasisMat(vertices, iQuad)
        cellM[:] += self.density * wt * numpy.dot(N.transpose(), N)
      feutils.assembleMat(self.M, cellM, cell, self.spaceDim)
    return


  def _calculateElasticityMat(self):
    """
    Calculate elasticity matrix.
    """
    return

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


  def _calculateBasisDerivMat(self):
    """
    Calculate matrix of derivatives of basis functions.
    """
    B = numpy.zeros( (self.spaceDim, self.spaceDim*self.numBasis),
                     dtype=numpy.float64)
    for iBasis in xrange(self.numBasis):
      for iDim in xrange(self.spaceDim):
        N[iDim, iBasis*self.spaceDim+iDim] = self.basis[iQuad, iBasis]
    return


# End of file 
