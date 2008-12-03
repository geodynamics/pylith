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
    self.quadrature = None
    self.gravityField = None
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


  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Do initialization.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    assert(None != self.cppHandle)
    self.cppHandle.normalizer = normalizer.cppHandle
    if None != self.gravityField:
      self.cppHandle.gravityField = self.gravityField.cppHandle
    
    self._logger.eventEnd(logEvent)
    return


  def timeStep(self, dt):
    """
    Set time step for advancing from time t to time t+dt.
    """
    assert(None != self.cppHandle)
    self.cppHandle.timeStep = dt
    return


  def stableTimeStep(self):
    """
    Get stable time step for advancing from time t to time t+dt.
    """
    logEvent = "%stimestep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    assert(None != self.cppHandle)
    dt = self.cppHandle.stableTimeStep

    self._logger.eventEnd(logEvent)
    return dt


  def useSolnIncr(self, flag):
    """
    Set behavior for using total field solution or incremental field solution.
    """
    logEvent = "%ssolnIncr" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    assert(None != self.cppHandle)
    self.cppHandle.useSolnIncr = flag

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


  def integrateResidual(self, residual, t, fields):
    """
    Integrate contributions to residual term at time t.
    """
    logEvent = "%sresidual" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    assert(None != self.cppHandle)
    self.cppHandle.integrateResidual(residual, t, fields.cppHandle,
                                     self.mesh.cppHandle,
				     self.mesh.coordsys.cppHandle)
    self._logger.eventEnd(logEvent)
    return


  def integrateJacobian(self, jacobian, t, fields):
    """
    Integrate contributions to Jacobian term at time t.
    """
    logEvent = "%sjacobian" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    assert(None != self.cppHandle)
    self.cppHandle.integrateJacobian(jacobian, t, fields.cppHandle,
                                     self.mesh.cppHandle)
    self._logger.eventEnd(logEvent)
    return


  def integrateResidualAssembled(self, residual, t, fields):
    """
    Integrate contributions to residual term at time t that do not
    require assembly over cells, vertices, or processors.
    """
    logEvent = "%sresidualAs" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    assert(None != self.cppHandle)
    self.cppHandle.integrateResidualAssembled(residual, t,
                                              fields.cppHandle,
                                              self.mesh.cppHandle,
                                              self.mesh.coordsys.cppHandle)
    self._logger.eventEnd(logEvent)
    return


  def integrateJacobianAssembled(self, jacobian, t, fields):
    """
    Integrate contributions to Jacobian term at time t that do not
    require assembly over cells, vertices, or processors.
    """
    logEvent = "%sjacobianAs" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    assert(None != self.cppHandle)
    self.cppHandle.integrateJacobianAssembled(jacobian, t,
                                              fields.cppHandle,
                                              self.mesh.cppHandle)
    self._logger.eventEnd(logEvent)
    return


  def poststep(self, t, dt, totalTime, fields):
    """
    Hook for doing stuff after advancing time step.
    """
    logEvent = "%spoststep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    assert(None != self.cppHandle)
    self.cppHandle.updateState(t, fields.cppHandle, self.mesh.cppHandle)

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
    logger.setClassName("FE Integrator")
    logger.initialize()

    events = ["verify",
              "init",
              "timestep",
              "solnIncr",
              "newJacobian",
              "residual",
              "jacobian",
              "residualAs",
              "jacobianAs",
              "poststep",
              "finalize"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._logger = logger
    return
  

# End of file 
