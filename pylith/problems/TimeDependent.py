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

## @file pylith/problems/TimeDependent.py
##
## @brief Python abstract base class for time dependent crustal
## dynamics problems.
##
## Factory: problem.

from Problem import Problem

# TimeDependent class
class TimeDependent(Problem):
  """
  Python abstract base class for time dependent crustal dynamics problems.

  Factory: problem.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Problem.Inventory):
    """
    Python object for managing TimeDependent facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing TimeDependent facilities and properties.
    ##
    ## \b Properties
    ## @li \b total_time Time duration for simulation.
    ## @li \b default_dt Default time step.
    ##
    ## \b Facilities
    ## @li \b formulation Formulation for solving PDE.
    ## @li \b checkpoint Checkpoint manager.

    import pyre.inventory

    from pyre.units.time import second
    totalTime = pyre.inventory.dimensional("total_time", default=0.0*second,
                          validator=pyre.inventory.greaterEqual(0.0*second))
    totalTime.meta['tip'] = "Time duration for simulation."

    dt = pyre.inventory.dimensional("default_dt", default=1.0*second,
                                 validator=pyre.inventory.greater(0.0*second))
    dt.meta['tip'] = "Default time step for simulation."

    from Implicit import Implicit
    formulation = pyre.inventory.facility("formulation",
                                          family="pde_formulation",
                                          factory=Implicit)
    formulation.meta['tip'] = "Formulation for solving PDE."

    from pylith.utils.CheckpointTimer import CheckpointTimer
    checkpointTimer = pyre.inventory.facility("checkpoint",
                                              family="checkpointer",
                                              factory=CheckpointTimer)
    checkpointTimer.meta['tip'] = "Checkpoint manager."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="timedependent"):
    """
    Constructor.
    """
    Problem.__init__(self, name)
    self._loggingPrefix = "PrTD "
    return


  def preinitialize(self, mesh):
    """
    Setup integrators for each element family (material/quadrature,
    bc/quadrature, etc.).
    """
    self._setupLogging()
    logEvent = "%spreinit" % self._loggingPrefix    
    self._logger.eventBegin(logEvent)
    
    self._info.log("Pre-initializing problem.")
    self.mesh = mesh
    self.formulation.preinitialize(mesh, self.materials, self.bc,
                                   self.interfaces)

    self._logger.eventEnd(logEvent)
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix    
    self._logger.eventBegin(logEvent)
    
    Problem.verifyConfiguration(self)
    self.formulation.verifyConfiguration()

    self._logger.eventEnd(logEvent)
    return
  

  def initialize(self):
    """
    Setup integrators for each element family (material/quadrature,
    bc/quadrature, etc.).
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self._info.log("Initializing problem.")
    self.formulation.initialize(self.dimension, self.totalTime, self.dt)

    self._logger.eventEnd(logEvent)
    return


  def run(self, app):
    """
    Solve time dependent problem.
    """
    logEvent = "%srun" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self._info.log("Solving problem.")
    self.checkpointTimer.toplevel = app # Set handle for saving state
    
    dt = self.formulation.stableTimeStep()
    if dt.value == 0.0: # If formulation returns 0.0, use default time step
      dt = self.dt
    t = self.formulation.startTime(self.dt)
    while t.value < self.totalTime.value:
      self._info.log("Main time loop, current time is t=%s" % t)
      
      # Checkpoint if necessary
      self.checkpointTimer.update(t)

      # Get stable time step
      dt = self.formulation.stableTimeStep()
      if dt.value == 0.0:
        # If formulation returns 0.0, use default time step
        dt = self.dt

      # Do stuff before advancing time step
      self._prestep(t, dt)

      # Advance in time
      self._step(t, dt)

      # Do stuff after advancing time step
      self._poststep(t, dt, self.totalTime)

      # Update time step
      t += dt

    self._logger.eventEnd(logEvent)
    return


  def finalize(self):
    """
    Cleanup after running problem.
    """
    logEvent = "%sfinalize" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self.formulation.finalize()

    self._logger.eventEnd(logEvent)
    return


  def checkpoint(self):
    """
    Save problem state for restart.
    """
    Problem.checkpoint()
    
    # Save state of this object
    raise NotImplementedError, "TimeDependent::checkpoint() not implemented."
  
    # Save state of children
    self.formulation.checkpoint()
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Problem._configure(self)
    self.totalTime = self.inventory.totalTime
    self.dt = self.inventory.dt
    self.formulation = self.inventory.formulation
    self.checkpointTimer = self.inventory.checkpointTimer
    return


  def _prestep(self, t, dt):
    """
    Hook for doing stuff before advancing time step.
    """
    self._info.log("Preparing to advance solution from time t=%s." % t)
    self.formulation.prestep(t, dt)
    return


  def _step(self, t, dt):
    """
    Advance to next time step.
    """
    self._info.log("Advancing solution from t=%s to t=%s." % (t, t+dt))    
    self.formulation.step(t, dt)
    return


  def _poststep(self, t, dt, totalTime):
    """
    Hook for doing stuff after advancing time step.
    """
    self._info.log("Finishing advancing solution to t=%s." % (t+dt))    
    self.formulation.poststep(t, dt, totalTime)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def problem():
  """
  Factory associated with TimeDependent.
  """
  return TimeDependent()


# End of file 
