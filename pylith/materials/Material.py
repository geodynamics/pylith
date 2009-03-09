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

## @file pylith/materials/Material.py
##

## @brief Python abstract base class for managing physical properties
## and state variables of a material.
##
## This implementation of a material associates both physical
## properties and a quadrature scheme with the material. Thus,
## applying different quadrature schemes within a region with the same
## physical property database requires two "materials", which can use
## the same database.
##
## Factory: material

from materials import ModuleMaterial
from pyre.components.Component import Component

# Material class
class Material(Component, ModuleMaterial):
  """
  Python material property manager.

  This implementation of a material associates both physical
  properties and a quadrature scheme with the material. Thus, applying
  different quadrature schemes within a region with the same physical
  property database requires two 'materials', which can use the same
  database.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing Material facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing Material facilities and properties.
    ##
    ## \b Properties
    ## @li \b id Material identifier (from mesh generator)
    ## @li \b name Name of material
    ## @li \b useInitialState Use initial state (true) or not (false).
    ##
    ## \b Facilities
    ## @li \b db Database of material property parameters
    ## @li \b quadrature Quadrature object for numerical integration
    ## @li \b dbInitialState Database for initial state.

    import pyre.inventory

    useInitialState = pyre.inventory.bool("use_initial_state", default=False)
    useInitialState.meta['tip'] = "Use initial state for material."

    id = pyre.inventory.int("id", default=0)
    id.meta['tip'] = "Material identifier (from mesh generator)."

    label = pyre.inventory.str("label", default="")
    label.meta['tip'] = "Name of material."

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    dbProperties = pyre.inventory.facility("properties_db",
                                           family="spatial_database",
                                           factory=SimpleDB)
    dbProperties.meta['tip'] = "Database for physical property parameters."

    dbInitialState = pyre.inventory.facility("initial_state_db",
                                           family="spatial_database",
                                           factory=SimpleDB)
    dbInitialState.meta['tip'] = "Database for initial state variables."

    from pylith.feassemble.quadrature.Quadrature import Quadrature
    quadrature = pyre.inventory.facility("quadrature", factory=Quadrature)
    quadrature.meta['tip'] = "Quadrature object for numerical integration."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="material"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="material")
    ModuleMaterial.__init__(self)
    self.output = None
    return


  def preinitialize(self):
    """
    Do pre-initialization setup.
    """
    self.quadrature.preinitialize()
    self._setupLogging()
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    if self.quadrature.spaceDim != self.dimension:
        raise ValueError, \
              "Quadrature scheme and material are incompatible.\n" \
              "Dimension for quadrature: %d\n" \
              "Dimension for material '%s': %d" % \
              (self.quadrature.spaceDim, self.label, self.dimension)
    
    self._logger.eventEnd(logEvent)
    return
  

  def getDataMesh(self):
    """
    Get mesh associated with data fields.
    """
    return (self.mesh, "material-id", self.id())


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Component._configure(self)
    self.id(self.inventory.id)
    self.label(self.inventory.label)
    self.dbProperties(self.inventory.dbProperties)
    if self.inventory.useInitialState:
      self.dbInitialState(self.inventory.dbInitialState)

    self.quadrature = self.inventory.quadrature
    return

  
  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    raise NotImplementedError("Please implement _createCppHandle() in " \
                              "derived class.")
  
  
  def _setupLogging(self):
    """
    Setup event logging.
    """
    if None == self._loggingPrefix:
      self._loggingPrefix = ""

    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.setClassName("FE Material")
    logger.initialize()

    events = ["verify",
              "init"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._logger = logger
    return
  

# End of file 
