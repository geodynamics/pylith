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

## @file pylith/problems/TimeStep.py
##
## @brief Python abstract base class for managing the time step size.
##
## Factory: time_step

from pyre.components.Component import Component

from pylith.utils.profiling import resourceUsageString

# TimeStep class
class TimeStep(Component):
  """
  Python abstract base class for managing the time step size.

  Factory: time_step.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="timestep"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="time_step")
    from pyre.units.time import second
    self.totalTime = 0.0*second
    self.dt = 0.0*second
    return


  def preinitialize(self):
    """
    Setup time step size algorithm.
    """
    self._setupLogging()
    logEvent = "%spreinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    self._logger.eventEnd(logEvent)
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self._logger.eventEnd(logEvent)
    return
  

  def initialize(self, normalizer):
    """
    Initialize time step algorithm.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    # Nondimensionalize time scales
    timeScale = normalizer.timeScale()
    self.totalTime = normalizer.nondimensionalize(self.totalTime, timeScale)
    self.dt = normalizer.nondimensionalize(self.dt, timeScale)

    self._logger.eventEnd(logEvent)
    return


  def numTimeSteps(self):
    """
    Get number of total time steps (or best guess if adaptive).
    """
    raise NotImplementedError("Please implement numTimeSteps().");
    return 0


  def timeStep(self, integrators):
    """
    Get stable time step for advancing forward in time.
    """
    # Default is to use the current time step.
    return self.dt
  

  def currentStep(self):
    """
    Get current time step size.
    """
    return self.dt
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    return


  def _setupLogging(self):
    """
    Setup event logging.
    """
    if not "_loggingPrefix" in dir(self):
      self._loggingPrefix = ""

    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.className("PDE TimeStep")
    logger.initialize()

    events = ["preinit",
              "verify",
              "init"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._logger = logger
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def time_step():
  """
  Factory associated with TimeStep.
  """
  return TimeStep()


# End of file 
