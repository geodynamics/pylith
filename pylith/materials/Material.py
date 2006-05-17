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
  """Python material property manager."""

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """Python object for managing Material facilities and properties."""
    
    ## @class Inventory
    ## Python object for managing Material facilities and properties.
    ##
    ## \b Properties
    ## @li \b id Material identifier (from mesh generator)
    ## @li \b name Name of material
    ##
    ## \b Facilities
    ## @li \b model Material model (constitutive equations)
    ## @li \b db Database of material property parameters

    import pyre.inventory

    id = pyre.inventory.int("id", default=0)
    id.meta['tip'] = "Material identifier (from mesh generator)."

    matname = pyre.inventory.str("name", default="")
    matname.meta['tip'] = "Name of material."

    from MaterialModel import MaterialModel
    model = pyre.inventory.facility("model", factory=MaterialModel,
                                    args=["model"])
    model.meta['tip'] = "Material model (constitutive equations)."

    from spatialdata.spatialdb.SpatialDB import SpatialDB
    db = pyre.inventory.facility("db", factory=SpatialDB,
                                 args=["db"])
    db.meta['tip'] = "Database of material properties."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def initialize(self):
    """Initialize material property manager."""

    self._info.log("Initializing material '%s'." % self.name)
    self.db.initialize()
    return


  def openDB(self):
    """Open material property database."""

    valNames = self.model.queryVals
    self._info.line("Material '%s' opening property database." % self.name)
    self._info.log("  Setting up query for values: %s." % valNames)
    self.db.open()
    self.db.queryVals(valNames)
    return


  def closeDB(self):
    """Close material property database."""

    self._info.log("Material '%s' closing property database." % self.name)
    self.db.close()
    return


  def __init__(self, name="material"):
    """Constructor."""
    
    Component.__init__(self, name, facility="material")
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """Setup members using inventory."""
    self.id = self.inventory.id
    self.matname = self.inventory.matname
    self.model = self.inventory.model
    self.db = self.inventory.db
    return

  
 # version
__id__ = "$Id$"

# End of file 
