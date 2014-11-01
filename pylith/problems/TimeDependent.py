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
# Copyright (c) 2010-2014 University of California, Davis
#
# See COPYING for license information.
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
    ## @li \b elastic_prestep Include a static calculation with elastic behavior before time stepping.
    ##
    ## \b Facilities
    ## @li \b formulation Formulation for solving PDE.
    ## @li \b progress_monitor Simple progress monitor via text file.
    ## @li \b checkpoint Checkpoint manager.

    import pyre.inventory

    elasticPrestep = pyre.inventory.bool("elastic_prestep", default=True)
    elasticPrestep.meta['tip'] = "Include a static calculation with elastic behavior before time stepping."

    from Implicit import Implicit
    formulation = pyre.inventory.facility("formulation", family="pde_formulation", factory=Implicit)
    formulation.meta['tip'] = "Formulation for solving PDE."

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
    Problem.__init__(self, name)
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
    self.formulation.preinitialize(mesh, self.materials, self.bc, self.interfaces, self.gravityField)
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    Problem.verifyConfiguration(self)
    self.formulation.verifyConfiguration()
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
    self.formulation.initialize(self.dimension, self.normalizer)
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
    
    # Elastic prestep
    if self.elasticPrestep:
      if 0 == comm.rank:
        self._info.log("Preparing for prestep with elastic behavior.")
      self._eventLogger.stagePush("Prestep")

      t = self.formulation.getStartTime()
      dt = self.formulation.getTimeStep()
      t -= dt

      # Limit material behavior to elastic regime
      for material in self.materials.components():
        material.useElasticBehavior(True)

      self.formulation.prestepElastic(t, dt)
      self._eventLogger.stagePop()

      if 0 == comm.rank:
        self._info.log("Computing prestep with elastic behavior.")
      self._eventLogger.stagePush("Step")
      self.formulation.step(t, dt)
      self._eventLogger.stagePop()

      if 0 == comm.rank:
        self._info.log("Finishing prestep with elastic behavior.")
      self._eventLogger.stagePush("Poststep")
      self.formulation.poststep(t, dt)
      self._eventLogger.stagePop()


    # Allow inelastic behavior
    for material in self.materials.components():
      material.useElasticBehavior(False)


    if (self.formulation.getTotalTime() > self.formulation.getStartTime()):
      self.progressMonitor.open()

    # Normal time loop
    t = self.formulation.getStartTime()
    timeScale = self.normalizer.timeScale()
    while t < self.formulation.getTotalTime():
      tsec = self.normalizer.dimensionalize(t, timeScale)
      tStart = self.normalizer.dimensionalize(self.formulation.getStartTime(), timeScale)
      tEnd = self.normalizer.dimensionalize(self.formulation.getTotalTime(), timeScale)
      self.progressMonitor.update(tsec, tStart, tEnd)

      self._eventLogger.stagePush("Prestep")
      if 0 == comm.rank:
        self._info.log("Main time loop, current time is t=%s" % tsec)
      
      # Checkpoint if necessary
      self.checkpointTimer.update(t)

      # Get time step for advancing in time
      dt = self.formulation.getTimeStep()
      dtsec = self.normalizer.dimensionalize(dt, timeScale)

      if 0 == comm.rank:
        self._info.log("Preparing to advance solution from time t=%s to t=%s." %\
                         (tsec, tsec+dtsec))
      self.formulation.prestep(t, dt)
      self._eventLogger.stagePop()

      if 0 == comm.rank:
        self._info.log("Advancing solution from t=%s to t=%s." % \
                         (tsec, tsec+dtsec))
      self._eventLogger.stagePush("Step")
      self.formulation.step(t, dt)
      self._eventLogger.stagePop()

      if 0 == comm.rank:
        self._info.log("Finishing advancing solution from t=%s to t=%s." % \
                         (tsec, tsec+dtsec))
      self._eventLogger.stagePush("Poststep")
      self.formulation.poststep(t, dt)
      self._eventLogger.stagePop()

      # Update time
      t += dt

    self.progressMonitor.close()
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
    self.elasticPrestep = self.inventory.elasticPrestep
    self.formulation = self.inventory.formulation
    self.progressMonitor = self.inventory.progressMonitor
    self.checkpointTimer = self.inventory.checkpointTimer
    return


# FACTORIES ////////////////////////////////////////////////////////////

def problem():
  """
  Factory associated with TimeDependent.
  """
  return TimeDependent()


# End of file 
