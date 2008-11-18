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
## of a material.
##
## This implementation of a material associates both physical
## properties and a quadrature scheme with the material. Thus,
## applying different quadrature schemes within a region with the same
## physical property database requires two "materials", which can use
## the same database.
##
## Factory: material

from pyre.components.Component import Component

# Material class
class Material(Component):
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
    ## @li \b initialStateDB Database for initial state.

    import pyre.inventory

    useInitialState = pyre.inventory.bool("use_initial_state", default=False)
    useInitialState.meta['tip'] = "Use initial state for material."

    id = pyre.inventory.int("id", default=0)
    id.meta['tip'] = "Material identifier (from mesh generator)."

    label = pyre.inventory.str("label", default="")
    label.meta['tip'] = "Name of material."

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    db = pyre.inventory.facility("db", family="spatial_database",
                                 factory=SimpleDB)
    db.meta['tip'] = "Database of material property parameters."

    initialStateDB = pyre.inventory.facility("initial_state_db",
                                              family="spatial_database",
                                              factory=SimpleDB)
    initialStateDB.meta['tip'] = "Database used for initial state."
    
    from pylith.feassemble.quadrature.Quadrature import Quadrature
    quadrature = pyre.inventory.facility("quadrature", factory=Quadrature)
    quadrature.meta['tip'] = "Quadrature object for numerical integration."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="material"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="material")
    self.cppHandle = None
    self.dimension = None
    self.output = None
    return


  def preinitialize(self):
    """
    Do pre-initialization setup.
    """
    self._createCppHandle()
    self.cppHandle.id = self.id
    self.cppHandle.label = self.label
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
  

  def initialize(self, mesh, totalTime, numTimeSteps):
    """
    Initialize material property manager.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self._info.log("Initializing material '%s'." % self.label)
    self.mesh = mesh
    assert(None != self.cppHandle)
    self.db.initialize()
    self.cppHandle.db = self.db.cppHandle
    if self.initialStateDB != None:
      self._info.log("Initializing initial state database.")
      self.initialStateDB.initialize()
      self.cppHandle.initialStateDB = self.initialStateDB.cppHandle
    self.cppHandle.initialize(mesh.cppHandle, mesh.coordsys.cppHandle,
                              self.quadrature.cppHandle)

    self._logger.eventEnd(logEvent)
    return


  def getDataMesh(self):
    """
    Get mesh associated with data fields.
    """
    return (self.mesh, "material-id", self.id)


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Component._configure(self)
    self.id = self.inventory.id
    self.label = self.inventory.label
    self.db = self.inventory.db
    self.quadrature = self.inventory.quadrature
    if self.inventory.useInitialState:
      self.initialStateDB = self.inventory.initialStateDB
    else:
      self.initialStateDB = None
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
