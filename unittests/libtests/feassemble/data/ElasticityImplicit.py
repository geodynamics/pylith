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

## @file unittests/libtests/feassemble/data/ElasticityImplicit.py

## @brief Python application for generating C++ data files for testing
## C++ ElasticityImplicit object.

from pyre.components.Component import Component

import numpy

# ----------------------------------------------------------------------

# ElasticityImplicit class
class ElasticityImplicit(Component):
  """
  Python application for generating C++ data files for testing C++
  ElasticityImplicit object.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="elasticityimplicit"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="formulation")
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def calculateResidual(self, integrator):
    """
    Calculate contribution to residual of operator for integrator.

    {r} = -[K]{u(t)}
    """
    K = integrator._calculateStiffnessMat()    

    residual = -numpy.dot(K, integrator.fieldT+integrator.fieldTIncr)
    return residual.flatten()


  def calculateJacobian(self, integrator):
    """
    Calculate contribution to Jacobian matrix of operator for integrator.

    [A] = [K]
    """
    K = integrator._calculateStiffnessMat()    

    jacobian = K
    return jacobian


# FACTORY //////////////////////////////////////////////////////////////
def formulation():
  return ElasticityImplicit()


# End of file 
