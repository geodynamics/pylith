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

## @file pylith/bc/BoundaryCondition.py
##
## @brief Python abstract base class for managing a boundary condition.
##
## This implementation of a boundary condition applies to a single
## face of an domain and associates both a quadrature scheme with a
## physical boundary condition. Thus, applying different quadrature
## schemes along a face with the same physical boundary condition
## requires two "bc", which can use the same database.
##
## Factory: boundary_condition

from pyre.components.Component import Component

# BoundaryCondition class
class BoundaryCondition(Component):
  """
  Python abstract base class for managing a boundary condition.

  This implementation of a boundary condition applies to a single
  face of an domain.

  Factory: boundary_condition
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing BoundaryCondition facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing BoundaryCondition facilities and properties.
    ##
    ## \b Properties
    ## @li \b id Object face identifier (from mesh generator)
    ## @li \b name Name identifier for object face
    ##
    ## \b Facilities
    ## @li \b db Database of boundary condition parameters

    import pyre.inventory

    id = pyre.inventory.int("id", default=0)
    id.meta['tip'] = "Object face identifier (from mesh generator)."

    label = pyre.inventory.str("label", default="")
    label.meta['tip'] = "Name identifier for object face."

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    db = pyre.inventory.facility("db", factory=SimpleDB,
                                 args=["db"])
    db.meta['tip'] = "Database of boundary condition parameters."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="boundarycondition"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="boundary_condition")
    self.cppHandle = None
    return


  def initialize(self, mesh):
    """
    Initialize boundary condition.
    """
    assert(None != self.cppHandle)
    self.db.initialize()
    self.cppHandle.id = self.id
    self.cppHandle.label = self.label
    self.cppHandle.db = self.db.cppHandle    
    self.mesh = mesh
    self.cppHandle.initialize(mesh.cppHandle, mesh.coordsys.cppHandle)
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
    return

  
# End of file 
