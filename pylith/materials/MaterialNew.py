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
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/materials/MaterialNew.py
##
## @brief Python abstract base class for managing physical properties
## and state variables of a material.
##
## Factory: material

from pylith.utils.PetscComponent import PetscComponent

# Validator for label
def validateLabel(value):
  """
  Validate descriptive label.
  """
  if 0 == len(value):
    raise ValueError("Descriptive label for material not specified.")
  return value


# MaterialNew class
class MaterialNew(PetscComponent):
  """
  Python material property manager.

  Factory: material
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing MaterialNew facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Material facilities and properties.
    ##
    ## \b Properties
    ## @li \b id Material identifier (from mesh generator)
    ## @li \b label Descriptive label for material.
    ##
    ## \b Facilities
    ## @li \b auxiliary_fields Discretization of auxiliary fields associated with material.
    ## @li \b db_auxiliary_fields Database for auxiliary fields associated with material.

    import pyre.inventory

    id = pyre.inventory.int("id", default=0)
    id.meta['tip'] = "Material identifier (from mesh generator)."

    label = pyre.inventory.str("label", default="", validator=validateLabel)
    label.meta['tip'] = "Descriptive label for material."

    from pylith.topology.AuxSubfield import subfieldFactory
    from pylith.utils.EmptyBin import EmptyBin
    auxFields = pyre.inventory.facilityArray("auxiliary_fields", itemFactory=subfieldFactory, factory=EmptyBin)
    auxFields.meta['tip'] = "Discretization of physical properties and state variables."

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    auxFieldsDB = pyre.inventory.facility("db_auxiliary_fields", family="spatial_database", factory=SimpleDB)
    auxFieldsDB.meta['tip'] = "Database for physical property parameters."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="material"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="material")
    self._createModuleObj()
    self.output = None
    return


  def preinitialize(self, mesh):
    """
    Do pre-initialization setup.
    """
    self._setupLogging()
    import weakref
    self.mesh = weakref.ref(mesh)
    from pylith.topology.topology import MeshOps_numMaterialCells
    self.ncells = MeshOps_numMaterialCells(mesh, self.id())
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    return


  def finalize(self):
    """
    Cleanup.
    """
    if not self.output is None:
      self.output.finalize()
    return


  def getDataMesh(self):
    """
    Get mesh associated with data fields.
    """
    return (self.mesh(), "material-id", self.id())


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    try:
      PetscComponent._configure(self)
      self.id(self.inventory.id)
      self.label(self.inventory.label)
      self.auxFieldsDB(self.inventory.auxFieldsDB)
      from pylith.utils.NullComponent import NullComponent
      if not isinstance(self.inventory.dbReferenceState, NullComponent):
        self.dbReferenceState(self.inventory.dbReferenceState)

    except ValueError, err:
      aliases = ", ".join(self.aliases)
      raise ValueError("Error while configuring material "
                       "(%s):\n%s" % (aliases, err.message))
    return


  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    raise NotImplementedError, \
          "Please implement _createModuleOb() in derived class."


  def _setupLogging(self):
    """
    Setup event logging.
    """
    if None == self._loggingPrefix:
      self._loggingPrefix = ""

    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.className("FE Material")
    logger.initialize()

    events = ["verify",
              "init"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._eventLogger = logger
    return


# End of file
