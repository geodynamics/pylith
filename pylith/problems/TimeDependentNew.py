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

## @file pylith/problems/TimeDependentNew.py
##
## @brief Python class for time dependent crustal
## dynamics problems.
##
## Factory: problem.

from ProblemNew import ProblemNew
from problems import TimeDependent as ModuleTimeDependent

# TimeDependentNew class
class TimeDependentNew(ProblemNew, ModuleTimeDependent):
  """
  Python class for time dependent crustal dynamics problems.

  Factory: problem.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(ProblemNew.Inventory):
    """
    Python object for managing TimeDependentNew facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing TimeDependentNew facilities and properties.
    ##
    ## \b Properties
    ## @li \b elastic_prestep Include a static calculation with elastic behavior before time stepping.
    ##
    ## \b Facilities
    ## @li \b initializer Problem initializer. 
    ## @li \b progress_monitor Simple progress monitor via text file.
    ## @li \b checkpoint Checkpoint manager.

    import pyre.inventory
    from pyre.units.time import year

    dtInitial = pyre.inventory.dimensional("initial_dt", default=1.0*year, validator=pyre.inventory.greater(0.0*year))
    dtInitial.meta['tip'] = "Initial time step."

    startTime = pyre.inventory.dimensional("start_time", default=0.0*year)
    startTime.meta['tip'] = "Start time for problem."

    totalTime = pyre.inventory.dimensional("total_time", default=0.0*year, validator=pyre.inventory.greaterEqual(0.0*year))
    totalTime.meta['tip'] = "Time duration of problem."

    initializer = pyre.inventory.facility("setup", family="initializer", factory=NullComponent)
    initializer.meta['tip'] = "Problem initializer."

    from ProgressMonitorTime import ProgressMonitorTime
    progressMonitor = pyre.inventory.facility("progress_monitor", family="progress_monitor", factory=ProgressMonitorTime)
    formulation.meta['tip'] = "Simple progress monitor via text file."

    from pylith.utils.CheckpointTimer import CheckpointTimer
    checkpointTimer = pyre.inventory.facility("checkpoint", family="checkpointer", factory=CheckpointTimer)
    checkpointTimer.meta['tip'] = "Checkpoint manager."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="timedependent"):
    """
    Constructor.
    """
    ProblemNew.__init__(self, name)
    ModuleTimeDependent.__init__(self)
    self._loggingPrefix = "PrTD "
    return


  def preinitialize(self, mesh):
    """
    Setup integrators for each element family (material/quadrature,
    bc/quadrature, etc.).
    """
    self._setupLogging()
    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()
    
    if 0 == comm.rank:
      self._info.log("Pre-initializing problem.")
    import weakref
    self.mesh = weakref.ref(mesh)
    ProblemNew.preinitialize(self, mesh)
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    ProblemNew.verifyConfiguration(self)
    return
  

  def initialize(self):
    """
    Setup integrators for each element family (material/quadrature,
    bc/quadrature, etc.).
    """
    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    if 0 == comm.rank:
      self._info.log("Initializing problem.")
    self.checkpointTimer.initialize(self.normalizer)
    ProblemNew.initialize(self)
    return


  def run(self, app):
    """
    Solve time dependent problem.
    """
    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    if 0 == comm.rank:
      self._info.log("Solving problem.")
    self.checkpointTimer.toplevel = app # Set handle for saving state
    
    # Pre-problem initialization
    if self.initializer:
      if 0 == comm.rank:
        self._info.log("")
      self._eventLogger.stagePush("Initialize")

    ModuleTimeDependent.create(self)
    ModuleTimeDependent.initialize()
    ModuleTimeDependent.solve()
    return


  def finalize(self):
    """
    Cleanup after running problem.
    """
    self.formulation.finalize()
    return


  def checkpoint(self):
    """
    Save problem state for restart.
    """
    ProblemNew.checkpoint()
    
    # Save state of this object
    raise NotImplementedError, "TimeDependentNew::checkpoint() not implemented."
  
    # Save state of children
    self.formulation.checkpoint()
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    ProblemNew._configure(self)
    ModuleTimeDependent.solver(self.inventory.solver)
    ModuleTimeDependent.initialTime(self.inventory.startTime)
    ModuleTimeDependent.timestep(self.inventory.timeStep)
    ModuleTimeDependent.duration(self.inventory.totalTime)

    self.initializer = self.inventory.initializer
    return


# FACTORIES ////////////////////////////////////////////////////////////

def problem():
  """
  Factory associated with TimeDependentNew.
  """
  return TimeDependentNew()


# End of file 
