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

## @brief Python material property manager.

from pyre.components.Component import Component

# Material class
class Material(Component):
  """
  Python material property manager.
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
    ##
    ## \b Facilities
    ## @li \b db Database of material property parameters

    import pyre.inventory

    id = pyre.inventory.int("id", default=0)
    id.meta['tip'] = "Material identifier (from mesh generator)."

    matname = pyre.inventory.str("name", default="")
    matname.meta['tip'] = "Name of material."

    from spatialdata.spatialdb.SpatialDB import SpatialDB
    db = pyre.inventory.facility("db", factory=SpatialDB,
                                 args=["db"])
    db.meta['tip'] = "Database of material property parameters."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def initialize(self):
    """
    Initialize material property manager.
    """
    self._info.log("Initializing material '%s'." % self.matname)
    self.db.initialize()
    return


  def __init__(self, name="material"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="material")
    self.cppHandle = None
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """Setup members using inventory."""
    self.id = self.inventory.id
    self.matname = self.inventory.matname
    self.db = self.inventory.db
    return

  
# End of file 
