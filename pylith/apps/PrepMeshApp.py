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

## @file pylith/apps/PrepMeshApp.py
##
## @brief Python pre-processing application to adjust topology of a
## mesh and partition it before running PyLith. This alleviates the
## need to redo these steps for every simulation when reusing the same
## faults.
##
## WARNING: This preprocessing application MUST be rerun if you change
## which faults are used in a problem.

from PetscApplication import PetscApplication

def faultFactory(name):
  """
  Factory for fault items.
  """
  from pyre.inventory import facility
  from pylith.faults.FaultCohesiveKin import FaultCohesiveKin
  return facility(name, family="fault", factory=FaultCohesiveKin)


# PrepMeshApp class
class PrepMeshApp(PetscApplication):
  """
  Python PrepMeshApp application.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  ## \b Properties
  ## @li None
  ##
  ## \b Facilities
  ## @li \b mesher Generates or imports the computational mesh.
  ## @li \b writer Writer for processed mesh.
  ## @li \b interfaces Interior surfaces with constraints or
  ##   constitutive models.

  import pyre.inventory
  from pylith.utils.EmptyBin import EmptyBin

  from pylith.topology.MeshImporter import MeshImporter
  mesher = pyre.inventory.facility("mesh_generator", family="mesh_generator",
                                   factory=MeshImporter)
  mesher.meta['tip'] = "Generates or imports the computational mesh."

  from pylith.meshio.MeshIOSieve import MeshIOSieve
  writer = pyre.inventory.facility("writer", family="mesh_io",
                                   factory=MeshIOSieve)
  writer.meta['tip'] = "Writer for processed mesh."

  interfaces = pyre.inventory.facilityArray("interfaces",
                                            itemFactory=faultFactory,
                                            factory=EmptyBin)
  interfaces.meta['tip'] = "Interior surfaces with constraints or " \
                           "constitutive models."

  from pylith.perf.MemoryLogger import MemoryLogger
  perfLogger = pyre.inventory.facility("perf_logger", family="perf_logger",
                                       factory=MemoryLogger)
  perfLogger.meta['tip'] = "Performance and memory logging."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="prepmeshapp"):
    """
    Constructor.
    """
    PetscApplication.__init__(self, name)
    self._loggingPrefix = "PrepMeshApp "
    return


  def main(self, *args, **kwds):
    """
    Run the application.
    """
    from pylith.utils.profiling import resourceUsageString
    
    self._debug.log(resourceUsageString())

    self._setupLogging()

    # Create mesh (adjust to account for interfaces (faults) if necessary)
    self._eventLogger.stagePush("Meshing")
    interfaces = None
    if "interfaces" in dir(self.problem):
      interfaces = self.problem.interfaces.components()
    mesh = self.mesher.create(self.problem.normalizer, interfaces)
    del interfaces
    del self.mesher
    self._debug.log(resourceUsageString())
    self._eventLogger.stagePop()

    self._eventLogger.stagePush("Output")
    writer.write(mesh)
    self._eventLogger.stagePop()


    # Cleanup
    self.perfLogger.logMesh('Mesh', mesh)
    self.compilePerformanceLog()
    if self.perfLogger.verbose:
      self.perfLogger.show()

    del mesh
    del self.problem
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

    self._eventLogger = logger
    return
  

# End of file 
