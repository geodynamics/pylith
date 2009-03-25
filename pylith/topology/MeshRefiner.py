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

## @file pylith/topology/MeshRefiner.py
##
## @brief Python manager for refining mesh in parallel.
##
## Factory: mesh_refiner.

from pyre.components.Component import Component

# MeshRefiner class
class MeshRefiner(Component):
  """
  Python manager for refining mesh in parallel.

  Factory: mesh_refiner
  """

  class Inventory(Component.Inventory):
    """
    Python object for managing RefineUniform facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing RefineUniform facilities and properties.
    ##
    ## \b Properties
    ## @li \b debug Write partition information to file.
    ##
    ## \b Facilities
    ## @li \b writer Data writer for for partition information.

    import pyre.inventory

    debug = pyre.inventory.bool("debug", default=False)
    debug.meta['tip'] = "Write partition information to file."

    from pylith.meshio.DataWriterVTK import DataWriterVTK
    dataWriter = pyre.inventory.facility("data_writer", factory=DataWriterVTK,
                                         family="output_data_writer")
    dataWriter.meta['tip'] = "Data writer for partition information."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="refiner"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="refiner")
    self.cppHandle = None
    return


  def refine(self, mesh):
    """
    Refine mesh.
    """
    return mesh


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.debug = self.inventory.debug
    self.dataWriter = self.inventory.dataWriter
    return


  def _createCppHandle(self):
    """
    Create handle to C++ object.
    """
    raise NotImplementedError("Please implement _createCppHandle().");
    return
  

  def _setupLogging(self):
    """
    Setup event logging.
    """
    self._loggingPrefix = "Refin "
    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.className("FE Refinement")
    logger.initialize()
    events = ["refine"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._logger = logger
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def mesh_refiner():
  """
  Factory associated with MeshRefiner.
  """
  return MeshRefiner()


# End of file 
