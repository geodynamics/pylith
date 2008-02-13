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

## @file pylith/feassemble/Constraint.py
##
## @brief Python abstract base class for constraints on operator
## actions with finite-elements.
##
## Factory: fe_constraint.

def implementsConstraint(obj):
  """
  Check whether object implements a constraint.
  """
  result = True
  available = dir(obj)
  required = ["preinitialize",
              "verifyConfiguration",
              "initialize",
              "timeStep",
              "setConstraintSizes",
              "setConstraints",
              "useSolnIncr",
              "setField",
              "poststep",
              "finalize"]
  for attr in required:
    if not attr in available:
      result = False
  return result



# Constraint class
class Constraint(object):
  """
  Python abstract base class for constraints on operator
  actions with finite-elements.

  Factory: constraint.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    self.cppHandle = None
    self.mesh = None
    return


  def preinitialize(self, mesh):
    """
    Setup constraint.
    """
    self.mesh = mesh
    self._setupLogging()
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
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


  def setConstraintSizes(self, field):
    """
    Set constraint sizes in field.
    """
    logEvent = "%ssetSizes" % self._loggingPrefix

    assert(None != self.cppHandle)
    self._logger.eventBegin(logEvent)
    self.cppHandle.setConstraintSizes(field, self.mesh.cppHandle)
    self._logger.eventEnd(logEvent)
    return


  def setConstraints(self, field):
    """
    Set constraints for field.
    """
    logEvent = "%sconstraints" % self._loggingPrefix

    assert(None != self.cppHandle)
    self._logger.eventBegin(logEvent)
    self.cppHandle.setConstraints(field, self.mesh.cppHandle)
    self._logger.eventEnd(logEvent)
    return


  def useSolnIncr(self, flag):
    """
    Set behavior for using total field solution or incremental field solution.
    """
    assert(None != self.cppHandle)
    self.cppHandle.useSolnIncr = flag
    return
  

  def setField(self, t, field):
    """
    Set constrained values in field at time t.
    """
    logEvent = "%ssetField" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    assert(None != self.cppHandle)
    self.cppHandle.setField(t.value, field, self.mesh.cppHandle)

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
    if None == self._loggingPrefix:
      self._loggingPrefix = ""

    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.setClassName("FE Constraint")
    logger.initialize()

    events = ["verify",
              "init",
              "setSizes",
              "constraints",
              "setField",
              "poststep",
              "finalize"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._logger = logger
    return
  

# End of file 
