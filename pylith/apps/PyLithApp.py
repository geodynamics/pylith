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

from mpi import Application

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

    from pylith.topology.MeshImporter import MeshImporter
    mesher = pyre.inventory.facility("mesh_generator", family="mesh_generator",
                                     factory=MeshImporter)
    mesher.meta['tip'] = "Generates or imports the computational mesh."

    from pylith.problems.TimeDependent import TimeDependent
    problem = pyre.inventory.facility("problem", family="problem",
                                      factory=TimeDependent)
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
    self._loggingPrefix = "PyLith "
    return


  def main(self, *args, **kwds):
    """
    Run the application.
    """
    from pylith.utils.profiling import resourceUsageString
    
    self.petsc.initialize()
    self._debug.log(resourceUsageString())

    self._setupLogging()
    self._logger.eventBegin("PyLith main")

    # Create mesh (adjust to account for interfaces (faults) if necessary)
    self._logger.stagePush("Meshing")
    interfaces = None
    if "interfaces" in dir(self.problem):
      interfaces = self.problem.interfaces.components()
    mesh = self.mesher.create(self.problem.normalizer, interfaces)
    self._debug.log(resourceUsageString())
    self._logger.stagePop()

    # Setup problem, verify configuration, and then initialize
    self._logger.stagePush("Setup")
    self.problem.preinitialize(mesh)
    self._debug.log(resourceUsageString())

    self.problem.verifyConfiguration()

    self.problem.initialize()
    self._debug.log(resourceUsageString())

    self._logger.stagePop()

    # Run problem
    self._logger.stagePush("Run")
    self.problem.run(self)
    self._debug.log(resourceUsageString())
    self._logger.stagePop()

    # Cleanup
    self._logger.stagePush("Finalize")
    self.problem.finalize()
    self._logger.stagePop()
    
    self._logger.eventEnd("PyLith main")
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

    import journal
    self._debug = journal.debug(self.name)
    return

  def _setupLogging(self):
    """
    Setup event logging.
    """
    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.className("PyLith")
    logger.initialize()
    logger.registerEvent("PyLith main")

    stages = ["Meshing",
              "Setup",
              "Run",
              "Finalize"]
    for stage in stages:
      logger.registerStage(stage)

    self._logger = logger
    return
  

# End of file 
