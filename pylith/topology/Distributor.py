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

## @file pylith/topology/Distributor.py
##
## @brief Python manager for distributing mesh among processors.
##
## Factory: mesh_distributor.

from pylith.utils.PetscComponent import PetscComponent
from topology import Distributor as ModuleDistributor

# Distributor class
class Distributor(PetscComponent, ModuleDistributor):
  """
  Python manager for distributing mesh among processors.

  Inventory

  \b Properties
  @li \b writePartition Write partition information to file.
  
  \b Facilities
  @li \b writer Data writer for for partition information.

  Factory: mesh_distributor
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory
    
  writePartition = pyre.inventory.bool("write_partition", default=False)
  writePartition.meta['tip'] = "Write partition information to file."
  
  from pylith.meshio.DataWriterVTK import DataWriterVTK
  dataWriter = pyre.inventory.facility("data_writer", factory=DataWriterVTK,
                                       family="data_writer")
  dataWriter.meta['tip'] = "Data writer for partition information."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="partitioner"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="partitioner")
    ModuleDistributor.__init__(self)
    return


  def distribute(self, mesh, normalizer):
    """
    Distribute a Mesh
    """
    self._setupLogging()
    logEvent = "%sdistribute" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    from pylith.topology.Mesh import Mesh
    newMesh = Mesh(mesh.dimension())
    ModuleDistributor.distribute(newMesh, mesh)

    #from pylith.utils.petsc import MemoryLogger
    #memoryLogger = MemoryLogger.singleton()

    #memoryLogger.stagePush(mesh.memLoggingStage)
    mesh.cleanup()
    #memoryLogger.stagePop()

    if self.writePartition:
      self.dataWriter.initialize(normalizer)
      ModuleDistributor.write(self.dataWriter, newMesh)

    self._eventLogger.eventEnd(logEvent)
    return newMesh


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)
    self.writePartition = self.inventory.writePartition
    self.dataWriter = self.inventory.dataWriter
    return


  def _setupLogging(self):
    """
    Setup event logging.
    """
    self._loggingPrefix = "Dist "
    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.className("FE Distribution")
    logger.initialize()
    events = ["distribute"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._eventLogger = logger
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def mesh_distributor():
  """
  Factory associated with Distributor.
  """
  return Distributor()


# End of file 
