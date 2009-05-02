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
  @li \b partitioner Name of mesh partitioner {"parmetis", "chaco"}.
  @li \b debug Write partition information to file.
  
  \b Facilities
  @li \b writer Data writer for for partition information.

  Factory: mesh_distributor
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory
    
  partitioner = pyre.inventory.str("partitioner", default="chaco",
                                   validator=pyre.inventory.choice(["chaco",
                                                                    "parmetis"]))
  partitioner.meta['tip'] = "Name of mesh partitioner."
  
  debug = pyre.inventory.bool("debug", default=False)
  debug.meta['tip'] = "Write partition information to file."
  
  from pylith.meshio.DataWriterVTKMesh import DataWriterVTKMesh
  dataWriter = pyre.inventory.facility("data_writer", factory=DataWriterVTKMesh,
                                       family="output_data_writer")
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
    self._logger.eventBegin(logEvent)

    from pylith.topology.Mesh import Mesh
    newMesh = Mesh(mesh.dimension())
    ModuleDistributor.distribute(newMesh, mesh, self.partitioner)

    if self.debug:
      self.dataWriter.initialize(normalizer)
      ModuleDistributor.write(self.dataWriter, newMesh)

    self._logger.eventEnd(logEvent)
    return newMesh


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)
    self.partitioner = self.inventory.partitioner
    self.debug = self.inventory.debug
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

    self._logger = logger
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def mesh_distributor():
  """
  Factory associated with Distributor.
  """
  return Distributor()


# End of file 
