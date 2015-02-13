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

## @file pylith/feassemble/Integrator.py
##
## @brief Python abstract base class for integration of operator
## actions with finite-elements.
##
## Factory: fe_integrator.

def implementsIntegrator(obj):
  """
  Check whether object implements an integrator.
  """
  result = True
  available = dir(obj)
  required = ["timeStep",
              "stableTimeStep",
              "integrateResidual",
              "integrateJacobian",
              "preinitialize",
              "verifyConfiguration",
              "initialize",
              "poststep",
              "finalize"]
  
  for attr in required:
    if not attr in available:
      result = False
  return result


# Integrator class
class Integrator(object):
  """
  Python abstract base class for integration of actions with
  finite-elements.

  Factory: integrator.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    self.mesh = None
    return


  def preinitialize(self, mesh):
    """
    Do pre-initialization setup.
    """
    self._setupLogging()
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    return


  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Do initialization.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    self.normalizer(normalizer)
    
    self._eventLogger.eventEnd(logEvent)
    return


  def poststep(self, t, dt, fields):
    """
    Hook for doing stuff after advancing time step.
    """
    logEvent = "%spoststep" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    self.updateStateVars(t, fields)

    self._eventLogger.eventEnd(logEvent)
    return


  def writeData(self, t, fields):
    """
    Hook for writing data at time t.
    """
    logEvent = "%swrite" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    self._eventLogger.eventEnd(logEvent)
    return


  def finalize(self):
    """
    Cleanup after time stepping.
    """
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _setupLogging(self):
    """
    Setup event logging.
    """
    if not "_loggingPrefix" in dir(self):
      self._loggingPrefix = ""

    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.className("FE Integrator")
    logger.initialize()

    events = ["preinit",
              "verify",
              "init",
              "poststep",
              "write",
              "finalize"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._eventLogger = logger
    return
  

# End of file 
