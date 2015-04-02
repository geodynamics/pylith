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

## @file pylith/apps/PetscApplication.py
##
## @brief Python PETSc application for creating an MPI application
## that uses PETSc.

from mpi import Application

# PetscApplication class
class PetscApplication(Application):
  """
  Python PETSc application for creating an MPI application that uses
  PETSc.

  Inventory:
    petsc Manager for PETSc options
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  # Dummy facility for passing options to PETSc
  from pylith.utils.PetscManager import PetscManager
  petsc = pyre.inventory.facility("petsc", family="petsc_manager", factory=PetscManager)
  petsc.meta['tip'] = "Manager for PETSc options."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="petscapp"):
    """
    Constructor.
    """
    Application.__init__(self, name)
    return


  def onComputeNodes(self, *args, **kwds):
    """
    Run the application in parallel on the compute nodes.
    """
    self.petsc.initialize()

    try:

      self.main(*args, **kwds)

    except Exception, err:
      self.cleanup() # Attempt to clean up memory.
      print "Fatal error. Calling MPI_Abort() to abort PyLith application."
      # Print stacktrace
      from traceback import print_exc
      print_exc()
      from pylith.mpi import mpi
      errorCode = -1
      mpi.mpi_abort(mpi.petsc_comm_world(), errorCode)

    self.cleanup()
    self.petsc.finalize()
    return
  

  def compilePerformanceLog(self):
    """
    Compile performance and memory information.
    """
    from pylith.utils.PetscComponent import PetscComponent
    for component in self.components():
      if isinstance(component, PetscComponent):
        component.compilePerformanceLog(self.perfLogger)
    return


  def cleanup(self):
    """
    Deallocate data structures.
    """
    from pylith.utils.PetscComponent import PetscComponent
    for component in self.components():
      if isinstance(component, PetscComponent):
        component.cleanup()
    self._cleanup()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Application._configure(self)
    self.petsc = self.inventory.petsc
    return


  def _cleanup(self):
    """
    Deallocate locally managed data structures.
    """
    return    


# End of file 
