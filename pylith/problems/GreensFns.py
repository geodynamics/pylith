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
# Copyright (c) 2010-2012 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/problems/GreensFns.py
##
## @brief Python abstract base class for time dependent crustal
## dynamics problems.
##
## Factory: problem.

from Problem import Problem

# GreensFns class
class GreensFns(Problem):
  """
  Python abstract base class for time dependent crustal dynamics problems.

  Factory: problem.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Problem.Inventory):
    """
    Python object for managing GreensFns facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing GreensFns facilities and properties.
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

  def __init__(self, name="greensfns"):
    """
    Constructor.
    """
    Problem.__init__(self, name)
    self._loggingPrefix = "PrGF "
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
    self.mesh = mesh
    self.formulation.preinitialize(mesh, self.materials, self.bc,
                                   self.interfaces, self.gravityField)
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
    Compute Green's functions associated with fault slip.
    """
    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    if 0 == comm.rank:
      self._info.log("Computing Green's functions.")
    self.checkpointTimer.toplevel = app # Set handle for saving state

    nimpulses = 0
    ipulse = 0;
    while ipulse < nimpulses:
      self._eventLogger.stagePush("Prestep")
      if 0 == comm.rank:
        self._info.log("Main loop, impulse %d of %d." % (ipulse+1, nimpulses))
      
      # Checkpoint if necessary
      self.checkpointTimer.update(t)

      if 0 == comm.rank:
        self._info.log("Preparing impulse %d of %d." % \
                         (ipulse+1, nimpulses))
      self.formulation.prestep(t, dt)
      self._eventLogger.stagePop()

      if 0 == comm.rank:
        self._info.log("Computing response to impulse %d of %d." %
                         (ipulse+1, nimpulses))
      self._eventLogger.stagePush("Step")
      self.formulation.step(t, dt)
      self._eventLogger.stagePop()

      if 0 == comm.rank:
        self._info.log("Finishing impulse %d of %d." % \
                         (ipulse+1, nimpulses))
      self._eventLogger.stagePush("Poststep")
      self.formulation.poststep(t, dt)
      self._eventLogger.stagePop()

      # Update time
      t += dt
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
    raise NotImplementedError, "GreensFns::checkpoint() not implemented."
  
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
  Factory associated with GreensFns.
  """
  return GreensFns()


# End of file 
