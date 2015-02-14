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

## @file unittests/libtests/feassemble/data/ElasticityApp.py

## @brief Python application for generating C++ data files for testing
## C++ elasticity integrator objects.

from IntegratorApp import IntegratorApp

import numpy
import feutils

# ----------------------------------------------------------------------

# ElasticityApp class
class ElasticityApp(IntegratorApp):
  """
  Python application for generating C++ data files for testing C++
  elasticity integrator objects.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(IntegratorApp.Inventory):
    """Python object for managing IntegratorApp facilities and properties."""

    ## @class Inventory
    ## Python object for managing ElasticityIntegrator facilities and
    ## properties.
    ##
    ## \b Properties
    ## @li \b useGravity Include gravitational body forces in residual.
    ##
    ## \b Facilities
    ## @li \b formulation Elasticity formulation.

    import pyre.inventory

    useGravity = pyre.inventory.bool("use_gravity", default=False)
    useGravity.meta['tip'] = "Include gravitational body forces in residual."

    from ElasticityImplicit import ElasticityImplicit
    formulation = pyre.inventory.facility("formulation",
                                          factory=ElasticityImplicit)
    formulation.meta['tip'] = "Elasticity formulation."


  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="elasticityapp"):
    """
    Constructor.
    """
    IntegratorApp.__init__(self, name)

    self.density = None
    self.lameMu = None
    self.lameLambda = None

    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members using inventory.
    """
    IntegratorApp._configure(self)
    self.useGravity = self.inventory.useGravity
    self.formulation = self.inventory.formulation
    return


  def _calculateResidual(self):
    """
    Calculate contribution to residual of operator for integrator.
    """
    self.valsResidual = self.formulation.calculateResidual(self)
    if self.useGravity:
      self.valsResidual += self._calculateGravityVec()
    return


  def _calculateJacobian(self):
    """
    Calculate contribution to Jacobian matrix of operator for integrator.
    """
    self.valsJacobian = self.formulation.calculateJacobian(self)
    return


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
        B = self._calculateBasisDerivMat(basisDeriv, iQuad)
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
      (jacobian, jacobianInv, jacobianDet, basisDeriv) = \
                 feutils.calculateJacobian(self.quadrature, vertices)
      for iQuad in xrange(self.numQuadPts):
        wt = self.quadWts[iQuad] * jacobianDet[iQuad]
        N = self._calculateBasisMat(iQuad)
        cellM[:] += self.density * wt * numpy.dot(N.transpose(), N)
      feutils.assembleMat(M, cellM, cell, self.spaceDim)
    return M


  def _calculateGravityVec(self):
    """
    Calculate body force vector.
    """
    gravityGlobal = numpy.zeros((self.numVertices*self.spaceDim),
                                dtype=numpy.float64)
    for cell in self.cells:
      gravityCell = numpy.zeros((self.spaceDim*self.numBasis))
      vertices = self.vertices[cell, :]
      (jacobian, jacobianInv, jacobianDet, basisDeriv) = \
                 feutils.calculateJacobian(self.quadrature, vertices)
      for iQuad in xrange(self.numQuadPts):
        wt = self.quadWts[iQuad] * jacobianDet[iQuad] * self.density
        for iBasis in xrange(self.numBasis):
          valI = wt * self.basis[iQuad, iBasis]
          for iDim in xrange(self.spaceDim):
            gravityCell[iDim + iBasis * self.spaceDim] += \
                             valI * self.gravityVec[iDim]
      feutils.assembleVec(gravityGlobal, gravityCell, cell, self.spaceDim)
    return gravityGlobal
    

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
      C2211 = self.lameLambda
      C2222 = lambda2mu
      C2212 = 0.0
      C1211 = 0.0
      C1222 = 0.0
      C1212 = 2.0*self.lameMu
      D = numpy.array([ [C1111, C1122, 0.5*C1112],
                        [C2211, C2222, 0.5*C2212],
                        [0.5*C1211, 0.5*C1222, 0.5*C1212] ],
                      dtype=numpy.float64)
    elif 3 == self.cellDim:
      lambda2mu = self.lameLambda + 2.0*self.lameMu
      C1111 = lambda2mu
      C1122 = self.lameLambda
      C1133 = self.lameLambda
      C1112 = 0.0
      C1123 = 0.0
      C1113 = 0.0
      C2211 = self.lameLambda
      C2222 = lambda2mu
      C2233 = self.lameLambda
      C2212 = 0.0
      C2223 = 0.0
      C2213 = 0.0
      C3311 = self.lameLambda
      C3322 = self.lameLambda
      C3333 = lambda2mu
      C3312 = 0.0
      C3323 = 0.0
      C3313 = 0.0
      C1211 = 0.0
      C1222 = 0.0
      C1233 = 0.0
      C1212 = 2.0*self.lameMu
      C1223 = 0.0
      C1213 = 0.0
      C2311 = 0.0
      C2322 = 0.0
      C2333 = 0.0
      C2312 = 0.0
      C2323 = 2.0*self.lameMu
      C2313 = 0.0
      C1311 = 0.0
      C1322 = 0.0
      C1333 = 0.0
      C1312 = 0.0
      C1323 = 0.0
      C1313 = 2.0*self.lameMu
      D = numpy.array([ [C1111, C1122, C1133, 0.5*C1112, 0.5*C1123, 0.5*C1113],
                        [C2211, C2222, C2233, 0.5*C2212, 0.5*C2223, 0.5*C2213],
                        [C3311, C3322, C3333, 0.5*C3312, 0.5*C3323, 0.5*C3313],
                        [0.5*C1211, 0.5*C1222, 0.5*C1233, 0.5*C1212, 0.5*C1223, 0.5*C1213],
                        [0.5*C2311, 0.5*C2322, 0.5*C2333, 0.5*C2312, 0.5*C2323, 0.5*C2313],
                        [0.5*C1311, 0.5*C1322, 0.5*C1333, 0.5*C1312, 0.5*C1323, 0.5*C1313] ],
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


  def _calculateBasisDerivMat(self, basisDeriv, iQuad):
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


# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = ElasticityApp()
  app.run()


# End of file 
