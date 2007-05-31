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

from IntegratorElasticity import IntegratorElasticity

import numpy

# ----------------------------------------------------------------------

# ElasticityExplicit class
class ElasticityExplicit(IntegratorElasticity):
  """
  Python application for generating C++ data files for testing C++
  ElasticityExplicit object.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="elasticityexplicit"):
    """
    Constructor.
    """
    IntegratorElasticity.__init__(self, name)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _calculateResidual(self):
    """
    Calculate contribution to residual of operator for integrator.

    {r} = (1/dt**2)[M](-{u(t+dt)} + 2 {u(t)} - {u(t-dt)}) -
          [K]{u(t)}
    """
    K = self._calculateStiffnessMat()    
    M = self._calculateMassMat()

    dispResult = -self.fieldTpdt + 2.0*self.fieldT - self.fieldTmdt
    self.valsResidual = 1.0/self.dt**2 * numpy.dot(M, dispResult) - \
                        numpy.dot(K, self.fieldT)
    return


  def _calculateJacobian(self):
    """
    Calculate contribution to Jacobian matrix of operator for integrator.

    [A] = (1/dt**2)[M]
    """
    M = self._calculateMassMat()

    self.valsJacobian = 1.0/self.dt**2 * M
    return


# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = ElasticityExplicit()
  app.run()


# End of file 
