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

## @file unittests/libtests/feassemble/data/ElasticityExplicit.py

## @brief Python application for generating C++ data files for testing
## C++ ElasticityExplicit object.

from pyre.components.Component import Component

import numpy
import feutils

# ----------------------------------------------------------------------

# ElasticityExplicit class
class ElasticityExplicit(Component):
  """
  Python application for generating C++ data files for testing C++
  ElasticityExplicit object.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="elasticityexplicit"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="formulation")
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def calculateResidual(self, integrator):
    """
    Calculate contribution to residual of operator for integrator.

    {r} = (1/dt**2)[M](-{u(t+dt)} + 2 {u(t)} - {u(t-dt)}) -
          [K]{u(t)}
    """
    M = integrator._calculateMassMat()
    Ml = self._lumpMatrix(M, integrator.numBasis, integrator.spaceDim)
    
    vel = (integrator.fieldT + integrator.fieldTIncr - integrator.fieldTmdt) / (2.0*integrator.dt)
    acc = (integrator.fieldTIncr - integrator.fieldT + integrator.fieldTmdt) / (integrator.dt**2) 
    dispAdj = integrator.fieldT + integrator.dt*integrator.normViscosity*vel
    
    residual = -Ml*acc
    residual.flatten()

    residual += self._elasticityResidual(integrator, dispAdj)
    return residual


  def calculateJacobian(self, integrator):
    """
    Calculate contribution to Jacobian matrix of operator for integrator.

    [A] = (1/dt**2)[M]
    """
    M = integrator._calculateMassMat()
    Ml = self._lumpMatrix(M, integrator.numBasis, integrator.spaceDim)

    jacobian = 1.0/integrator.dt**2 * Ml
    return jacobian.flatten()


  def _lumpMatrix(self, matrix, numBasis, spaceDim):
    """
    Lump matrix.
    """
    (nrows,ncols) = matrix.shape
    assert(numBasis * spaceDim == nrows)
    assert(nrows == ncols)
    vector = numpy.zeros( (ncols), dtype=numpy.float64)

    for iBasis in xrange(numBasis):
      for iDim in xrange(spaceDim):
        v = 0.0
        for jBasis in xrange(numBasis):
          v += matrix[iBasis*spaceDim+iDim,jBasis*spaceDim+iDim]
        vector[iBasis*spaceDim+iDim] = v
    return vector


  def _elasticityResidual(self, integrator, dispAdj):
    """
    Calculate action for elasticity.
    """
    residual = numpy.zeros( (integrator.numBasis*integrator.spaceDim),
                            dtype=numpy.float64)

    # Matrix of elasticity values
    D = integrator._calculateElasticityMat()
    
    for cell in integrator.cells:
      cellR = numpy.zeros( (integrator.spaceDim*integrator.numBasis, 1),
                           dtype=numpy.float64)
      vertices = integrator.vertices[cell, :]
      (jacobian, jacobianInv, jacobianDet, basisDeriv) = \
          feutils.calculateJacobian(integrator.quadrature, vertices)

      for iQuad in xrange(integrator.numQuadPts):
        wt = integrator.quadWts[iQuad] * jacobianDet[iQuad]
        BL0 = integrator._calculateBasisDerivMatLinear0(basisDeriv, iQuad)
        BL1 = integrator._calculateBasisDerivMatLinear1(basisDeriv, iQuad, dispAdj)
        BL = BL0 + BL1
        strain = integrator._calculateStrain(basisDeriv, iQuad, dispAdj)
        S = numpy.dot(D, strain.transpose())
        cellR -= wt * numpy.dot(BL.transpose(), S)
      
      feutils.assembleVec(residual, cellR.flatten(), cell, integrator.spaceDim)
    return residual


# FACTORY //////////////////////////////////////////////////////////////
def formulation():
  return ElasticityExplicit()


# End of file 
