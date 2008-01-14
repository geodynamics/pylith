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

from pyre.components.Component import Component

# Distributor class
class Distributor(Component):
  """
  Python manager for distributing mesh among processors.

  Factory: mesh_distributor
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing Distributor facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Distributor facilities and properties.
    ##
    ## \b Properties
    ## @li \b partitioner Name of mesh partitioner {"parmetis", "chaco"}.
    ## @li \b debug Write partition information to file.
    ##
    ## \b Facilities
    ## @li \b output Output manager for partition information

    import pyre.inventory

    partitioner = pyre.inventory.str("partitioner", default="chaco",
                                     validator=pyre.inventory.choice(["chaco",
                                                                      "parmetis"]))
    partitioner.meta['tip'] = "Name of mesh partitioner."

    debug = pyre.inventory.bool("debug", default=False)
    debug.meta['tip'] = "Write partition information to file."

    from pylith.meshio.OutputManager import OutputManager
    output = pyre.inventory.facility("output", factory=OutputManager,
                                     family="output_manager")
    output.meta['tip'] = "Output manager for partition information."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="partitioner"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="partitioner")
    self.cppHandle = None
    self.debug = False
    return


  def distribute(self, mesh):
    """
    Distribute a Mesh
    """
    self._setupLogging()
    logEvent = "%sdistribute" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self._createCppHandle()
    
    from Mesh import Mesh
    newMesh = Mesh()
    newMesh.cppHandle = self.cppHandle.distribute(mesh.cppHandle,
                                                  self.partitioner)
    newMesh.coordsys = mesh.coordsys

    if self.debug:
      self.output.initialize()
      self.cppHandle.write(self.output.cppHandle,
                           newMesh.cppHandle, newMesh.coordsys.cppHandle)

    self._logger.eventEnd(logEvent)
    return newMesh


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.partitioner = self.inventory.partitioner
    self.debug = self.inventory.debug
    self.output = self.inventory.output
    return


  def _createCppHandle(self):
    """
    Create handle to C++ object.
    """
    if None == self.cppHandle:
      import pylith.topology.topology as bindings
      self.cppHandle = bindings.Distributor()
    return
  

  def _setupLogging(self):
    """
    Setup event logging.
    """
    self._loggingPrefix = "Dist "
    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.setClassName("FE Distribution")
    logger.initialize()
    events = ["distribute"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._logger = logger
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def mesh_partitioner():
  """
  Factory associated with Distributor.
  """
  return Distributor()


# End of file 
