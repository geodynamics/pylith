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
from geometry import Scene as ModuleScene

# Scene class
class Scene(Component, ModuleScene):

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """Python object for managing Scene facilities and properties."""

    import pyre.inventory

    from Bar import Bar
    bar = pyre.inventory.facility("bar", family="shape", factory=Bar)
    bar.meta['tip'] = "Bar in scene."

    from Sphere import Sphere
    sphere = pyre.inventory.facility("sphere", family="shape", factory=Sphere)
    sphere.meta['tip'] = "Sphere in scene."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="scene"):
    """
    Constructor.
    """
    Component.__init__(self, name)
    ModuleScene.__init__(self)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Component._configure(self)

    self.bar(self.inventory.bar)
    self.sphere(self.inventory.sphere)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def scene():
  return Scene()


# End of file 
