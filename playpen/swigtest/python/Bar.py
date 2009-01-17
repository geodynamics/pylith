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

from Shape import Shape
from geometry import Bar as ModuleBar

# Bar class
class Bar(Shape, ModuleBar):

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Shape.Inventory):
    """Python object for managing Bar facilities and properties."""

    import pyre.inventory

    length = pyre.inventory.float("length", default=1.0)
    length.meta['tip'] = "Length of bar."

    width = pyre.inventory.float("width", default=1.0)
    width.meta['tip'] = "Width of bar."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="bar"):
    """
    Constructor.
    """
    Shape.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Shape._configure(self)

    self.length(self.inventory.length)
    self.width(self.inventory.width)
    return


  def _createModuleObj(self):
    ModuleBar.__init__(self)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def shape():
  return Bar()


# End of file 
