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

## @file pylith/feassemble/ElasticityExplicitLgDeform.py
##
## @brief Python object for explicit time integration of dynamic
## elasticity equation using finite-elements with large rigidy body
## motion and small strains.
##
## Factory: integrator

from IntegratorElasticityLgDeform import IntegratorElasticityLgDeform
from feassemble import ElasticityExplicitLgDeform as ModuleElasticityExplicitLgDeform

# ElasticityExplicitLgDeform class
class ElasticityExplicitLgDeform(IntegratorElasticityLgDeform,
                                 ModuleElasticityExplicitLgDeform):
  """
  Python object for explicit time integration of dynamic elasticity
  equation using finite-elements with large rigid body motions and
  small strains.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticityexplicit"):
    """
    Constructor.
    """
    IntegratorElasticityLgDeform.__init__(self, name)
    ModuleElasticityExplicitLgDeform.__init__(self)
    self._loggingPrefix = "ElEx "
    return


  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Do initialization.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    IntegratorElasticityLgDeform.initialize(self, totalTime, numTimeSteps, normalizer)
    ModuleElasticityExplicitLgDeform.initialize(self, self.mesh)
    self._initializeOutput(totalTime, numTimeSteps, normalizer)
    
    self._eventLogger.eventEnd(logEvent)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def integrator():
  """
  Factory associated with ElasticityExplicitLgDeform.
  """
  return ElasticityExplicitLgDeform()


# End of file 
