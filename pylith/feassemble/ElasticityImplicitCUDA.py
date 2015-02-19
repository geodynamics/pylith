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

## @file pylith/feassemble/ElasticityImplicitCUDA.py
##
## @brief Python object for implicit time integration of dynamic
## elasticity equation using finite-elements.
##
## Factory: integrator

from IntegratorElasticity import IntegratorElasticity
from feassemble import ElasticityImplicitCUDA as ModuleElasticityImplicit

# ElasticityImplicit class
class ElasticityImplicitCUDA(IntegratorElasticity, ModuleElasticityImplicit):
  """
  Python object for implicit time integration of elasticity
  equation using finite-elements.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticityimplicit"):
    """
    Constructor.
    """
    IntegratorElasticity.__init__(self, name)
    ModuleElasticityImplicit.__init__(self)
    self._loggingPrefix = "ElIm "
    return


  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Do initialization.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    IntegratorElasticity.initialize(self, totalTime, numTimeSteps, normalizer)
    ModuleElasticityImplicit.initialize(self, self.mesh())
    self._initializeOutput(totalTime, numTimeSteps, normalizer)
    
    self._eventLogger.eventEnd(logEvent)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def integrator():
  """
  Factory associated with ElasticityImplicitCUDA.
  """
  return ElasticityImplicitCUDA()


# End of file 
