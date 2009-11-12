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

## @file unittests/libtests/feassemble/data/ElasticityImplicitGrav.py

## @brief Python application for generating C++ data files for testing
## C++ ElasticityImplicitGrav object.

from ElasticityImplicit import ElasticityImplicit

import numpy
import feutils

# ----------------------------------------------------------------------

# ElasticityImplicitGrav class
class ElasticityImplicitGrav(ElasticityImplicit):
  """
  Python application for generating C++ data files for testing C++
  ElasticityImplicitGrav object.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="elasticityimplicitgrav"):
    """
    Constructor.
    """
    ElasticityImplicit.__init__(self, name)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def calculateResidual(self, integrator):
    """
    Calculate contribution to residual of operator for integrator.

    {r} = -[K]{u(t)}
    """
    residual = ElasticityImplicit.calculateResidual(self, integrator)
    gravityGlobal = self._calculateGravity(integrator)
    residual += gravityGlobal.reshape(residual.shape)
    return residual


  def _calculateGravity(self, integrator):
    """
    Calculate body force vector.
    """
    gravityGlobal = numpy.zeros(( integrator.numVertices*integrator.spaceDim ),
                                dtype=numpy.float64)
    for cell in integrator.cells:
      gravityCell = numpy.zeros(integrator.spaceDim*integrator.numBasis)
      vertices = integrator.vertices[cell, :]
      (jacobian, jacobianInv, jacobianDet, basisDeriv) = \
                 feutils.calculateJacobian(integrator.quadrature, vertices)
      for iQuad in xrange(integrator.numQuadPts):
        wt = integrator.quadWts[iQuad] * jacobianDet[iQuad] * integrator.density
        for iBasis in xrange(integrator.numBasis):
          valI = wt * integrator.basis[iQuad, iBasis]
          for iDim in xrange(integrator.spaceDim):
            gravityCell[iDim + iBasis * integrator.spaceDim] += \
                             valI * integrator.gravityVec[iDim]
      feutils.assembleVec(gravityGlobal, gravityCell, cell, integrator.spaceDim)
    return gravityGlobal
    

# FACTORY //////////////////////////////////////////////////////////////
def formulation():
  return ElasticityImplicitGrav()


# End of file 
