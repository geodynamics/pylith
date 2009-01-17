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

from pyre.components.Component import Component

# Shape class
class Shape(Component):

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """Python object for managing Shape facilities and properties."""

    import pyre.inventory

    color = pyre.inventory.str("color", default="black")
    color.meta['tip'] = "Color of object."

    id = pyre.inventory.int("id", default=1)
    id.meta['tip'] = "Id number for bar."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="shape"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="shape")
    self._createModuleObj()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Component._configure(self)

    self.color(self.inventory.color)
    self.id(self.inventory.id)
    return


  def _createModuleObj(self):
    raise NotImplementedError("_createModuleObj() not implemented.")
    return
  

# End of file 
