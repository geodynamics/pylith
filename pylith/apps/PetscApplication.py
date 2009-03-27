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
  petsc = pyre.inventory.facility("petsc", family="petsc_manager",
                                  factory=PetscManager)
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
    self.main(*args, **kwds)
    self.cleanup()
    self.petc.finalize()
    return
  

  def cleanup(self):
    """
    Deallocate data structures.
    """
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
