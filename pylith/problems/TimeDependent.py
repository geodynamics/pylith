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
    ## @li \b dimension Spatial dimension of problem space.
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

    dimension = pyre.inventory.int("dimension", default=3,
                                   validator=pyre.inventory.choice([1,2,3]))
    dimension.meta['tip'] = "Spatial dimension of problem space."

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
    return


  def initialize(self, mesh):
    """
    Setup integrators for each element family (material/quadrature,
    bc/quadrature, etc.).
    """
    self._info.log("Initializing problem.")

    if self.dimension != mesh.dimension():
        raise ValueError, \
              "Spatial dimension of problem is '%d' but mesh contains cells " \
              "for spatial dimension '%d'." % \
              (self.dimension, mesh.dimension)
    self.mesh = mesh
    self.formulation.initialize(mesh, self.materials, self.bc,
                                self.interfaces, self.dimension, self.dt)
    return


  def run(self, app):
    """
    Solve time dependent problem.
    """
    self._info.log("Solving problem.")
    self.checkpointTimer.toplevel = app # Set handle for saving state
    
    from pyre.units.time import second
    t = 0.0*second
    while t.value <= self.totalTime.value:
      self._info.log("Main time loop, t=%s" % t)
      
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
      self._step(dt)

      # Do stuff after advancing time step
      self._poststep(t+dt)

      # Update time step
      t += dt
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
    self.dimension = self.inventory.dimension
    self.formulation = self.inventory.formulation
    self.checkpointTimer = self.inventory.checkpointTimer
    return


  def _prestep(self, t, dt):
    """
    Hook for doing stuff before advancing time step.
    """
    self.formulation.prestep(t, dt)
    return


  def _step(self, dt):
    """
    Advance to next time step.
    """
    self.formulation.step(dt)
    return


  def _poststep(self, t):
    """
    Hook for doing stuff after advancing time step.
    """
    self.formulation.poststep(t)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def problem():
  """
  Factory associated with TimeDependent.
  """
  return TimeDependent()


# End of file 
