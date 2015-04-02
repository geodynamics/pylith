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

    {r} = -Sum(wt * [BL]^T [S})
    """
    import feutils

    residual = numpy.zeros( (integrator.spaceDim*integrator.numVertices),
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
        cellR -= wt * numpy.dot(BL.transpose(), S)
      
      feutils.assembleVec(residual, cellR.flatten(), cell, integrator.spaceDim)

    return residual


# FACTORY //////////////////////////////////////////////////////////////
def formulation():
  return ElasticityImplicitLgDeform()


# End of file 
