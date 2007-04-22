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

## @file pylith/PyLithApp.py
##
## @brief Python PyLith application

#from mpi.Application import Application
from pyre.applications.Script import Script as Application

# PyLithApp class
class PyLithApp(Application):
  """
  Python PyLithApp application.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Application.Inventory):
    """
    Python object for managing PyLithApp facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing PyLithApp facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b mesher Generates or imports the computational mesh.
    ## @li \b problem Computational problem to solve
    ## @li \b petsc Manager for PETSc options

    import pyre.inventory

    from pyre.units.time import second
    totalTime = pyre.inventory.dimensional("total_time", default=0.0*second,
                          validator=pyre.inventory.greaterEqual(0.0*second))
    totalTime.meta['tip'] = "Time duration for simulation."

    from pylith.topology.MeshImporter import MeshImporter
    mesher = pyre.inventory.facility("mesh_generator", family="mesh_generator",
                                     factory=MeshImporter)
    mesher.meta['tip'] = "Generates or imports the computational mesh."

    from pylith.problems.EqDeformation import EqDeformation
    problem = pyre.inventory.facility("problem", family="problem",
                                      factory=EqDeformation)
    problem.meta['tip'] = "Computational problem to solve."

    # Dummy facility for passing options to PETSc
    from pylith.utils.PetscManager import PetscManager
    petsc = pyre.inventory.facility("petsc", family="petsc_manager",
                                    factory=PetscManager)
    petsc.meta['tip'] = "Manager for PETSc options."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="pylithapp"):
    """
    Constructor.
    """
    Application.__init__(self, name)
    return


  def main(self):
    """
    Run the application.
    """
    self.petsc.initialize()

    # Create mesh (adjust to account for faults if necessary)
    faults = None
    if faults in dir(self.problem):
      faults = self.problem.faults
    mesh = self.mesher.create(faults)

    # Initialize problem and then run
    self.problem.initialize(mesh)
    self.problem.run(self)
    
    self.petsc.finalize()
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Application._configure(self)
    self.mesher = self.inventory.mesher
    self.problem = self.inventory.problem
    self.petsc = self.inventory.petsc
    return


# End of file 
