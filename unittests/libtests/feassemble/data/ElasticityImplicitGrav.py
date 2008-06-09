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

from IntegratorElasticity import IntegratorElasticity

import numpy
import feutils
import pdb

# ----------------------------------------------------------------------

# ElasticityImplicitGrav class
class ElasticityImplicitGrav(IntegratorElasticity):
  """
  Python application for generating C++ data files for testing C++
  ElasticityImplicitGrav object.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="elasticityimplicitgrav"):
    """
    Constructor.
    """
    pdb.set_trace()
    IntegratorElasticity.__init__(self, name)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _calculateResidual(self):
    """
    Calculate contribution to residual of operator for integrator.

    {r} = -[K]{u(t)}
    """
    pdb.set_trace()
    K = self._calculateStiffnessMat()    
    gravityGlobal = self._calculateGravity()

    self.valsResidual = -numpy.dot(K, self.fieldTpdt) + gravityGlobal
    return


  def _calculateJacobian(self):
    """
    Calculate contribution to Jacobian matrix of operator for integrator.

    [A] = [K]
    """
    K = self._calculateStiffnessMat()    

    self.valsJacobian = K
    return


  def _calculateGravity(self):
    """
    Calculate body force vector.
    """
    pdb.set_trace()
    gravityGlobal = numpy.zeros(self.spaceDim*self.numVertices,
                                dtype=numpy.float64)
    for cell in self.cells:
      gravityCell = numpy.zeros(self.spaceDim*self.numBasis)
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
    


# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = ElasticityImplicitGrav()
  app.run()


# End of file 
