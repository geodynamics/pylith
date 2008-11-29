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
    ## None
    ##
    ## \b Facilities
    ## @li \b formulation Formulation for solving PDE.
    ## @li \b checkpoint Checkpoint manager.

    import pyre.inventory

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
                                   self.interfaces, self.gravityField)

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
    self.normalizer.initialize()
    self.formulation.initialize(self.dimension, self.normalizer)

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
    
    t = self.formulation.getStartTime()
    while t < self.formulation.getTotalTime():
      self._info.log("Main time loop, current time is t=%s" % t)
      
      # Checkpoint if necessary
      self.checkpointTimer.update(t)

      # Get time step for advancing in time
      dt = self.formulation.getTimeStep()

      self._info.log("Preparing to advance solution from time t=%s to t=%s." % (t, t+dt))
      self.formulation.prestep(t, dt)

      self._info.log("Advancing solution from t=%s to t=%s." % (t, t+dt))
      self.formulation.step(t, dt)

      self._info.log("Finishing advancing solution from t=%s to t=%s." % (t, t+dt))
      self.formulation.poststep(t, dt)

      # Update time
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
    self.formulation = self.inventory.formulation
    self.checkpointTimer = self.inventory.checkpointTimer
    return


# FACTORIES ////////////////////////////////////////////////////////////

def problem():
  """
  Factory associated with TimeDependent.
  """
  return TimeDependent()


# End of file 
