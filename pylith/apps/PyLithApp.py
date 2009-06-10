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

## @file pylith/apps/PyLithApp.py
##
## @brief Python PyLith application

from PetscApplication import PetscApplication

# PyLithApp class
class PyLithApp(PetscApplication):
  """
  Python PyLithApp application.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscApplication.Inventory):
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

    from pylith.perf.MemoryLogger import MemoryLogger
    perfLogger = pyre.inventory.facility("perf_logger", family="perf_logger",
                                         factory=MemoryLogger)
    perfLogger.meta['tip'] = "Performance and memory logging."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="pylithapp"):
    """
    Constructor.
    """
    PetscApplication.__init__(self, name)
    self._loggingPrefix = "PyLith "
    return


  def main(self, *args, **kwds):
    """
    Run the application.
    """
    from pylith.utils.profiling import resourceUsageString
    
    self._debug.log(resourceUsageString())

    self._setupLogging()

    # Create mesh (adjust to account for interfaces (faults) if necessary)
    self._logger.stagePush("Meshing")
    interfaces = None
    if "interfaces" in dir(self.problem):
      interfaces = self.problem.interfaces.components()
    mesh = self.mesher.create(self.problem.normalizer, interfaces)
    self.perfLogger.logMesh('Mesh', mesh)
    del interfaces
    del self.mesher
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
    self.problem.run(self)
    self._debug.log(resourceUsageString())

    # Cleanup
    self._logger.stagePush("Finalize")
    self.problem.finalize()
    self._logger.stagePop()

    del mesh
    del self.problem

    self.compilePerformanceLog()
    if self.perfLogger.verbose: self.perfLogger.show()
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    PetscApplication._configure(self)
    self.mesher = self.inventory.mesher
    self.problem = self.inventory.problem
    self.perfLogger = self.inventory.perfLogger

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

    self._logger = logger
    return
  

# End of file 
