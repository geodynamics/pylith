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
              "useSolnIncr",
              "integrateResidual",
              "integrateJacobian",
              "integrateResidualAssembled",
              "integrateJacobianAssembled",
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
    self.gravityField = None
    return


  def preinitialize(self, mesh):
    """
    Do pre-initialization setup.
    """
    self._setupLogging()
    return


  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Do initialization.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self.normalizer(normalizer)
    if None != self.gravityField:
      self.gravityField(self.gravityField)
    
    self._logger.eventEnd(logEvent)
    return


  def poststep(self, t, dt, totalTime, fields):
    """
    Hook for doing stuff after advancing time step.
    """
    logEvent = "%spoststep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self.updateState(t, fields)

    self._logger.eventEnd(logEvent)
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
              "init",
              "poststep",
              "finalize"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._logger = logger
    return
  

# End of file 
