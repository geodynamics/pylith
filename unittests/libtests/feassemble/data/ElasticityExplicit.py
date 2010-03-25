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
          [K]{u(t)}
    """
    K = integrator._calculateStiffnessMat()    
    M = integrator._calculateMassMat()

    vel = (integrator.fieldT + integrator.fieldTIncr - integrator.fieldTmdt) / (2.0*integrator.dt)
    acc = (integrator.fieldTIncr - integrator.fieldT + integrator.fieldTmdt) / (integrator.dt**2) 
    residual = -numpy.dot(M, acc) - numpy.dot(K, integrator.fieldT)
    return residual.flatten()


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
