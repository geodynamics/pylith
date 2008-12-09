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
from geometry import Sphere as ModuleSphere

# Sphere class
class Sphere(Shape, ModuleSphere):

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Shape.Inventory):
    """Python object for managing Sphere facilities and properties."""

    import pyre.inventory

    radius = pyre.inventory.float("radius", default=1.0)
    radius.meta['tip'] = "Radius of sphere."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="sphere"):
    """
    Constructor.
    """
    Shape.__init__(self, name)
    ModuleSphere.__init__(self)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Shape._configure(self)

    self.radius(self.inventory.radius)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def shape():
  return Sphere()


# End of file 
