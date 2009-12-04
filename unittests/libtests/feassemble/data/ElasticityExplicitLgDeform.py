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

## @file unittests/libtests/feassemble/data/ElasticityExplicit.py

## @brief Python application for generating C++ data files for testing
## C++ ElasticityExplicit object.

from pyre.components.Component import Component

import numpy

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
          Sum(wt * [BL]^T [S])
    """
    import feutils

    # Calculate action for inertia
    M = integrator._calculateMassMat()
    dispResult = integrator.fieldT - integrator.fieldTmdt
    residual = 1.0/integrator.dt**2 * numpy.dot(M, dispResult)
    residual = residual.flatten()

    # Calculate action for elasticity

    # Matrix of elasticity values
    D = integrator._calculateElasticityMat()
    
    for cell in integrator.cells:
      cellR = numpy.zeros( (integrator.spaceDim*integrator.numBasis, 1),
                           dtype=numpy.float64)
      vertices = integrator.vertices[cell, :]
      (jacobian, jacobianInv, jacobianDet, basisDeriv) = \
          feutils.calculateJacobian(integrator.quadrature, vertices)
      fieldT = integrator.fieldT
      for iQuad in xrange(integrator.numQuadPts):
        wt = integrator.quadWts[iQuad] * jacobianDet[iQuad]
        BL0 = integrator._calculateBasisDerivMatLinear0(basisDeriv, iQuad)
        BL1 = integrator._calculateBasisDerivMatLinear1(basisDeriv, iQuad, fieldT)
        BL = BL0 + BL1
        strain = integrator._calculateStrain(basisDeriv, iQuad, fieldT)
        S = numpy.dot(D, strain.transpose())
        cellR -= wt * numpy.dot(BL.transpose(), S)
      
      feutils.assembleVec(residual, cellR.flatten(), cell, integrator.spaceDim)

    return residual


  def calculateJacobian(self, integrator):
    """
    Calculate contribution to Jacobian matrix of operator for integrator.

    [A] = (1/dt**2)[M]
    """
    M = integrator._calculateMassMat()

    jacobian = 1.0/integrator.dt**2 * M
    return jacobian


# FACTORY //////////////////////////////////////////////////////////////
def formulation():
  return ElasticityExplicit()


# End of file 
