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

## @file unittests/libtests/feassemble/data/ElasticityImplicit.py

## @brief Python application for generating C++ data files for testing
## C++ ElasticityImplicit object.

from IntegratorElasticity import IntegratorElasticity

import numpy

# ----------------------------------------------------------------------

# ElasticityImplicit class
class ElasticityImplicit(IntegratorElasticity):
  """
  Python application for generating C++ data files for testing C++
  ElasticityImplicit object.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="elasticityimplicit"):
    """
    Constructor.
    """
    IntegratorElasticity.__init__(self, name)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _calculateResidual(self):
    """
    Calculate contribution to residual of operator for integrator.

    {r} = -[K]{u(t)}
    """
    K = self._calculateStiffnessMat()    

    self.valsResidual = -numpy.dot(K, self.fieldTpdt)
    return


  def _calculateJacobian(self):
    """
    Calculate contribution to Jacobian matrix of operator for integrator.

    [A] = [K]
    """
    K = self._calculateStiffnessMat()    

    self.valsJacobian = K
    return


# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = ElasticityImplicit()
  app.run()


# End of file 
