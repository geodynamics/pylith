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

## @file pylith/problems/TimeStep.py
##
## @brief Python abstract base class for managing the time step size.
##
## Factory: time_step

from pylith.utils.PetscComponent import PetscComponent
from pylith.utils.profiling import resourceUsageString

# TimeStep class
class TimeStep(PetscComponent):
  """
  Python abstract base class for managing the time step size.

  Factory: time_step.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing TimeStepUniform facilities and properties.
    """

    ## @class Inventory
    ## Python abstract base class for managing TimeStep facilities and properties.
    ##
    ## \b Properties
    ## @li \b total_time Time duration for simulation.
    ## @li \b start_time Starting time for simulation.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    from pyre.units.time import second
    totalTime = pyre.inventory.dimensional("total_time", default=0.0*second,
                          validator=pyre.inventory.greaterEqual(0.0*second))
    totalTime.meta['tip'] = "Time duration for simulation."

    startTime = pyre.inventory.dimensional("start_time", default=0.0*second)
    startTime.meta['tip'] = "Time duration for simulation."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="timestep"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="time_step")
    from pyre.units.time import second
    self.timeScale = 1.0*second
    self.totalTime = 0.0*second
    self.startTime = 0.0*second
    self.dt = 0.0*second
    self.totalTimeN = 0.0 # Nondimensionalized total time
    self.startTimeN = 0.0 # Nondimensionalized start time
    self.dtN = 0.0 # Nondimenionalized time step
    return


  def preinitialize(self):
    """
    Setup time step size algorithm.
    """
    self._setupLogging()
    logEvent = "%spreinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)
    
    self._eventLogger.eventEnd(logEvent)
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    self._eventLogger.eventEnd(logEvent)
    return
  

  def initialize(self, normalizer):
    """
    Initialize time step algorithm.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    # Nondimensionalize time scales
    timeScale = normalizer.timeScale()
    self.totalTimeN = normalizer.nondimensionalize(self.totalTime, timeScale)
    self.startTimeN = normalizer.nondimensionalize(self.startTime, timeScale)
    self.dtN = normalizer.nondimensionalize(self.dt, timeScale)
    self.timeScale = timeScale

    self._eventLogger.eventEnd(logEvent)
    return


  def numTimeSteps(self):
    """
    Get number of total time steps (or best guess if adaptive).
    """
    raise NotImplementedError("Please implement numTimeSteps().");
    return 0


  def timeStep(self, mesh, integrators):
    """
    Get stable time step for advancing forward in time.
    """
    # Default is to use the current time step.
    return self.dtN
  

  def currentStep(self):
    """
    Get current time step size.
    """
    return self.dtN
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)
    self.totalTime = self.inventory.totalTime
    self.startTime = self.inventory.startTime
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

    self._eventLogger = logger
    return
  

  def _stableTimeStep(self, mesh, integrators):
    """
    Get stable time step.
    """
    dtStable = 1.0e+30
    for integrator in integrators:
      dt = integrator.stableTimeStep(mesh)
      if dt < dtStable:
        dtStable = dt
    import pylith.mpi.mpi as mpi
    comm = mesh.comm()
    dtStableAll = mpi.allreduce_scalar_double(dtStable, mpi.mpi_min(), comm.handle)
    return dtStableAll



# FACTORIES ////////////////////////////////////////////////////////////

def time_step():
  """
  Factory associated with TimeStep.
  """
  return TimeStep()


# End of file 
