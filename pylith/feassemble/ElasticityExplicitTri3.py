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

## @file pylith/feassemble/ElasticityExplicitTri3.py
##
## @brief Python object for explicit time integration of dynamic
## elasticity equation using finite-elements.
##
## Factory: integrator

from IntegratorElasticity import IntegratorElasticity
from feassemble import ElasticityExplicitTri3 as ModuleElasticityExplicitTri3

# ElasticityExplicitTri3 class
class ElasticityExplicitTri3(IntegratorElasticity, ModuleElasticityExplicitTri3):
  """
  Python object for explicit time integration of dynamic elasticity
  equation using finite-elements.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticityexplicittri3"):
    """
    Constructor.
    """
    IntegratorElasticity.__init__(self, name)
    ModuleElasticityExplicitTri3.__init__(self)
    self._loggingPrefix = "ElEx "
    return


  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Do initialization.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    IntegratorElasticity.initialize(self, totalTime, numTimeSteps, normalizer)
    ModuleElasticityExplicitTri3.initialize(self, self.mesh)
    self._initializeOutput(totalTime, numTimeSteps, normalizer)
    
    self._eventLogger.eventEnd(logEvent)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def integrator():
  """
  Factory associated with ElasticityExplicitTri3.
  """
  return ElasticityExplicitTri3()


# End of file 
