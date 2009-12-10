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

## @file pylith/feassemble/ElasticityImplicitLgDeform.py
##
## @brief Python object for implicit time integration of dynamic
## elasticity equation using finite-elements.
##
## Factory: integrator

from IntegratorElasticityLgDeform import IntegratorElasticityLgDeform
from feassemble import ElasticityImplicitLgDeform as ModuleElasticityImplicitLgDeform

# ElasticityImplicitLgDeform class
class ElasticityImplicitLgDeform(IntegratorElasticityLgDeform,
                                 ModuleElasticityImplicitLgDeform):
  """
  Python object for implicit time integration of elasticity
  equation using finite-elements.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticityimplicit"):
    """
    Constructor.
    """
    IntegratorElasticityLgDeform.__init__(self, name)
    ModuleElasticityImplicitLgDeform.__init__(self)
    self._loggingPrefix = "ElIm "
    return


  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Do initialization.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    IntegratorElasticityLgDeform.initialize(self, totalTime, numTimeSteps,
                                            normalizer)
    ModuleElasticityImplicitLgDeform.initialize(self, self.mesh)
    self._initializeOutput(totalTime, numTimeSteps, normalizer)
    
    self._eventLogger.eventEnd(logEvent)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def integrator():
  """
  Factory associated with ElasticityImplicitLgDeform.
  """
  return ElasticityImplicitLgDeform()


# End of file 
