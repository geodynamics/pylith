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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/friction/FrictionModel.py
##

## @brief Python abstract base class for managing physical properties
## and state variables of a fault constitutive model.
##
## This implementation of a fault constitutive model associates both
## physical properties and a quadrature scheme with a fault. Thus,
## applying different quadrature schemes within a region with the same
## physical property database requires two "friction models", which
## can use the same database.
##
## Factory: friction_model

from pylith.utils.PetscComponent import PetscComponent

# Validator for label
def validateLabel(value):
  """
  Validate descriptive label.
  """
  if 0 == len(value):
    raise ValueError("Descriptive label for friction model not specified.")
  return value


# FrictionModel class
class FrictionModel(PetscComponent):
  """
  Python friction model property manager.

  This implementation of a fault constitutive model associates both
  physical properties and a quadrature scheme with a fault. Thus,
  applying different quadrature schemes within a region with the same
  physical property database requires two 'friction models', which can
  use the same database.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing FrictionModel facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing FrictionModel facilities and properties.
    ##
    ## \b Properties
    ## @li \b name Name of friction model.
    ##
    ## \b Facilities
    ## @li \b db_properties Database of material property parameters
    ## @li \b db_initial_state Database for initial state.

    import pyre.inventory

    label = pyre.inventory.str("label", default="", validator=validateLabel)
    label.meta['tip'] = "Descriptive label for friction model."

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    dbProperties = pyre.inventory.facility("db_properties",
                                           family="spatial_database",
                                           factory=SimpleDB)
    dbProperties.meta['tip'] = "Database for physical property parameters."

    from pylith.utils.NullComponent import NullComponent
    dbInitialState = pyre.inventory.facility("db_initial_state",
                                           family="spatial_database",
                                           factory=NullComponent)
    dbInitialState.meta['tip'] = "Database for initial state variables."

    from pylith.perf.MemoryLogger import MemoryLogger
    perfLogger = pyre.inventory.facility("perf_logger", family="perf_logger",
                                         factory=MemoryLogger)
    perfLogger.meta['tip'] = "Performance and memory logging."

  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="frictionmodel"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="frictionmodel")
    self._createModuleObj()
    return


  def finalize(self):
    """
    Cleanup.
    """
    self._modelMemoryUse()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    try:
      PetscComponent._configure(self)
      self.label(self.inventory.label)
      self.dbProperties(self.inventory.dbProperties)
      from pylith.utils.NullComponent import NullComponent
      if not isinstance(self.inventory.dbInitialState, NullComponent):
        self.dbInitialState(self.inventory.dbInitialState)

      self.perfLogger = self.inventory.perfLogger
    except ValueError, err:
      aliases = ", ".join(self.aliases)
      raise ValueError("Error while configuring friction model "
                       "(%s):\n%s" % (aliases, err.message))
    return

  
  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    raise NotImplementedError, \
          "Please implement _createModuleOb() in derived class."


  def _modelMemoryUse(self):
    """
    Model allocated memory.
    """
    self.perfLogger.logFrictionModel('Friction', self)
    self.perfLogger.logField('Friction', self.propertiesField())
    self.perfLogger.logField('Friction', self.stateVarsField())
    return


  def _setupLogging(self):
    """
    Setup event logging.
    """
    if None == self._loggingPrefix:
      self._loggingPrefix = ""

    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.className("FE FrictionModel")
    logger.initialize()

    events = []
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._eventLogger = logger
    return
  

# End of file 
