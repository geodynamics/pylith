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


  def preinitialize(self):
    """
    Do pre-initialization setup.
    """
    self._createCppHandle()
    self.cppHandle.id = self.id
    self.cppHandle.label = self.label
    self.quadrature.preinitialize()
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    if self.quadrature.spaceDim != self.dimension:
        raise ValueError, \
              "Quadrature scheme and material are incompatible.\n" \
              "Dimension for quadrature: %d\n" \
              "Dimension for material '%s': %d" % \
              (self.quadrature.spaceDim, self.label, self.dimension)
    
    # :TODO: Make sure mesh contains material (need to account for the
    # fact that any given processor may only have a subset of the
    # materials)
    return
  

  def initialize(self, mesh):
    """
    Initialize material property manager.
    """
    self._info.log("Initializing material '%s'." % self.label)
    assert(None != self.cppHandle)
    self.db.initialize()
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

  
  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    raise NotImplementedError("Please implement _createCppHandle() in " \
                              "derived class.")
  
  
# End of file 
