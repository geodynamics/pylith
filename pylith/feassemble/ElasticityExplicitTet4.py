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

## @file pylith/feassemble/ElasticityExplicitTet4.py
##
## @brief Python object for explicit time integration of dynamic
## elasticity equation using finite-elements.
##
## Factory: integrator

from IntegratorElasticity import IntegratorElasticity
from feassemble import ElasticityExplicitTet4 as ModuleElasticityExplicitTet4

# ElasticityExplicitTet4 class
class ElasticityExplicitTet4(IntegratorElasticity, ModuleElasticityExplicitTet4):
  """
  Python object for explicit time integration of dynamic elasticity
  equation using finite-elements.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticityexplicittet4"):
    """
    Constructor.
    """
    IntegratorElasticity.__init__(self, name)
    ModuleElasticityExplicitTet4.__init__(self)
    self._loggingPrefix = "ElEx "
    return


  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Do initialization.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    IntegratorElasticity.initialize(self, totalTime, numTimeSteps, normalizer)
    ModuleElasticityExplicitTet4.initialize(self, self.mesh())
    self._initializeOutput(totalTime, numTimeSteps, normalizer)
    
    self._eventLogger.eventEnd(logEvent)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _verifyConfiguration(self):
    ModuleElasticityExplicitTet4.verifyConfiguration(self, self.mesh())
    return


# FACTORIES ////////////////////////////////////////////////////////////

def integrator():
  """
  Factory associated with ElasticityExplicitTet4.
  """
  return ElasticityExplicitTet4()


# End of file 
