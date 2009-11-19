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

## @file unittests/libtests/feassemble/data/ElasticityImplicitLgDeform.py

## @brief Python application for generating C++ data files for testing
## C++ ElasticityImplicitLgDeform object.

from ElasticityImplicit import ElasticityImplicit

import numpy

# ----------------------------------------------------------------------

# ElasticityImplicitLgDeform class
class ElasticityImplicitLgDeform(ElasticityImplicit):
  """
  Python application for generating C++ data files for testing C++
  ElasticityImplicitLgDeform object.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="elasticityimplicitlgdeform"):
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
    import feutils

    residual = numpy.zeros( (integrator.spaceDim*integrator.numVertices, 1),
                            dtype=numpy.float64)

    # Matrix of elasticity values
    D = integrator._calculateElasticityMat()
    
    for cell in integrator.cells:
      cellR = numpy.zeros( (integrator.spaceDim*integrator.numBasis, 1),
                           dtype=numpy.float64)
      vertices = integrator.vertices[cell, :]
      (jacobian, jacobianInv, jacobianDet, basisDeriv) = \
          feutils.calculateJacobian(integrator.quadrature, vertices)
      fieldTpdt = integrator.fieldT + integrator.fieldTIncr
      for iQuad in xrange(integrator.numQuadPts):
        wt = integrator.quadWts[iQuad] * jacobianDet[iQuad]
        BL0 = integrator._calculateBasisDerivMatLinear0(basisDeriv, iQuad)
        BL1 = integrator._calculateBasisDerivMatLinear1(basisDeriv, iQuad, fieldTpdt)
        BL = BL0 + BL1
        strain = integrator._calculateStrain(basisDeriv, iQuad, fieldTpdt)
        S = numpy.dot(D, strain.transpose())
        print "BL:",BL.shape
        print "S:",S.shape
        cellR -= wt * numpy.dot(BL.transpose(), S)
        print "S",S
        print cellR
      feutils.assembleVec(residual, cellR, cell, integrator.spaceDim)

    return residual


# FACTORY //////////////////////////////////////////////////////////////
def formulation():
  return ElasticityImplicitLgDeform()


# End of file 
