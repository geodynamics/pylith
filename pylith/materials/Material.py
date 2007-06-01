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
    ##
    ## \b Facilities
    ## @li \b db Database of material property parameters
    ## @li \b quadrature Quadrature object for numerical integration

    import pyre.inventory

    id = pyre.inventory.int("id", default=0)
    id.meta['tip'] = "Material identifier (from mesh generator)."

    label = pyre.inventory.str("label", default="")
    label.meta['tip'] = "Name of material."

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    db = pyre.inventory.facility("db", factory=SimpleDB,
                                 args=["db"])
    db.meta['tip'] = "Database of material property parameters."
    
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
    return


  def initialize(self, mesh):
    """
    Initialize material property manager.
    """
    self._info.log("Initializing material '%s'." % self.label)

    if self.dimension != self.quadrature.cell.cellDim:
      raise ValueError, \
            "Quadrature is incompatible with material.\n" \
            "Dimensions for quadrature: %d, dimensions for material: %d" % \
            (self.quadrature.cell.cellDim, self.dimension)
    if self.dimension != mesh.dimension():
      raise ValueError, \
            "Material is incompatible with mesh.\n" \
            "Dimensions for mesh: %d, dimensions for material: %d" % \
            (mesh.dimension(), self.dimension)

    self.db.initialize()
    assert(None != self.cppHandle)
    self.cppHandle.id = self.id
    self.cppHandle.label = self.label
    self.cppHandle.db = self.db.cppHandle
    self.cppHandle.initialize(mesh.cppHandle, mesh.coordsys.cppHandle,
                              self.quadrature.cppHandle)
    return


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
    return

  
# End of file 
