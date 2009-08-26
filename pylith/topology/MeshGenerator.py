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

## @file pylith/topology/MeshGenerator.py
##
## @brief Python abstract base class for mesh generator.
##
## Factory: mesh_generator.

from pylith.utils.PetscComponent import PetscComponent

# MeshGenerator class
class MeshGenerator(PetscComponent):
  """
  Python abstract base class for mesh generator.

  Factory: mesh_generator
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing MeshGenerator facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing MeshGenerator facilities and properties.
    ##
    ## \b Properties
    ## @li \b debug Debugging flag for mesh.
    ## @li \b interpolate Build intermediate mesh topology elements (if true)
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    debug = pyre.inventory.bool("debug", default=False)
    debug.meta['tip'] = "Debugging flag for mesh."

    interpolate = pyre.inventory.bool("interpolate", default=False)
    interpolate.meta['tip'] = "Build intermediate mesh topology elements"


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="meshgenerator"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="meshgenerator")
    self.debug = False
    self.interpolate = False
    return


  def create(self, normalizer, faults=None):
    """
    Generate a Mesh.
    """
    raise NotImplementedError("MeshGenerator.create() not implemented.")
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)
    self.debug = self.inventory.debug
    self.interpolate = self.inventory.interpolate
    return


  def _adjustTopology(self, mesh, interfaces):
    """
    Adjust topology for interface implementation.
    """
    logEvent = "%sadjTopo" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)
    
    #self._info.activate()
    #mesh.view("===== MESH BEFORE ADJUSTING TOPOLOGY =====")

    if not interfaces is None:
      firstFaultVertex = 0
      firstFaultCell   = 0
      for interface in interfaces:
        self._info.log("Counting vertices for fault '%s'." % interface.label())
        nvertices = interface.numVertices(mesh)
        firstFaultCell += nvertices
        if interface.useLagrangeConstraints():
          firstFaultCell += nvertices
      for interface in interfaces:
        nvertices = interface.numVertices(mesh)
        self._info.log("Adjusting topology for fault '%s' with %d vertices." % \
                         (interface.label(), nvertices))
        firstFaultVertex, firstFaultCell = \
            interface.adjustTopology(mesh, firstFaultVertex, firstFaultCell)
        
    #mesh.view("===== MESH AFTER ADJUSTING TOPOLOGY =====")
    #self._info.deactivate()

    self._eventLogger.eventEnd(logEvent)
    return
  

  def _setupLogging(self):
    """
    Setup event logging.
    """
    if not "_loggingPrefix" in dir(self):
      self._loggingPrefix = ""

    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.className("Mesh Generator")
    logger.initialize()

    events = ["create",
              "adjTopo"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._eventLogger = logger
    return
  

# End of file
