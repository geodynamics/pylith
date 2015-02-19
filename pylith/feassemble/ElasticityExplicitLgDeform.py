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
    ModuleElasticityExplicitLgDeform.initialize(self, self.mesh())
    self._initializeOutput(totalTime, numTimeSteps, normalizer)
    
    self._eventLogger.eventEnd(logEvent)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _verifyConfiguration(self):
    ModuleElasticityExplicitLgDeform.verifyConfiguration(self, self.mesh())
    return


# FACTORIES ////////////////////////////////////////////////////////////

def integrator():
  """
  Factory associated with ElasticityExplicitLgDeform.
  """
  return ElasticityExplicitLgDeform()


# End of file 
