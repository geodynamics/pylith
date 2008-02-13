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
              "updateState",
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
    self.quadrature = None
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
    logEvent = "%sverify" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    assert(None != self.cppHandle)
    self.cppHandle.verifyConfiguration(self.mesh.cppHandle)

    self._logger.eventEnd(logEvent)
    return


  def initialize(self, totalTime, numTimeSteps):
    """
    Do initialization.
    """
    return


  def timeStep(self, dt):
    """
    Set time step for advancing from time t to time t+dt.
    """
    assert(None != self.cppHandle)
    self.cppHandle.timeStep = dt.value
    return


  def stableTimeStep(self):
    """
    Get stable time step for advancing from time t to time t+dt.
    """
    logEvent = "%stimestep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    assert(None != self.cppHandle)
    from pyre.units.time import second
    dt = self.cppHandle.stableTimeStep*second

    self._logger.eventEnd(logEvent)
    return dt


  def useSolnIncr(self, flag):
    """
    Set behavior for using total field solution or incremental field solution.
    """
    assert(None != self.cppHandle)
    self.cppHandle.useSolnIncr = flag
    return
  

  def integrateResidual(self, residual, t, fields):
    """
    Integrate contributions to residual term at time t.
    """
    logEvent = "%sresidual" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    assert(None != self.cppHandle)
    self.cppHandle.integrateResidual(residual, t.value, fields.cppHandle,
                                     self.mesh.cppHandle)
    self._logger.eventEnd(logEvent)
    return


  def needNewJacobian(self):
    """
    Returns true if we need to recompute Jacobian matrix for operator,
    false otherwise.
    """
    logEvent = "%snewJacobian" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    assert(None != self.cppHandle)
    flag = self.cppHandle.needNewJacobian
    self._logger.eventEnd(logEvent)
    return flag


  def integrateJacobian(self, jacobian, t, fields):
    """
    Integrate contributions to Jacobian term at time t.
    """
    logEvent = "%sjacobian" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    assert(None != self.cppHandle)
    self.cppHandle.integrateJacobian(jacobian, t.value, fields.cppHandle,
                                     self.mesh.cppHandle)
    self._logger.eventEnd(logEvent)
    return


  def updateState(self, t, fields):
    """
    Update state variables as needed.
    """
    logEvent = "%sstate" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    assert(None != self.cppHandle)
    self.cppHandle.updateState(t.value, fields.cppHandle, self.mesh.cppHandle)

    self._logger.eventEnd(logEvent)
    return
    

  def poststep(self, t, dt, totalTime, fields):
    """
    Hook for doing stuff after advancing time step.
    """
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
    logger.setClassName("FE Integrator")
    logger.initialize()

    events = ["verify",
              "init",
              "timestep",
              "residual",
              "newJacobian",
              "jacobian",
              "state",
              "poststep",
              "finalize"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._logger = logger
    return
  

# End of file 
