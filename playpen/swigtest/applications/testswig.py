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

from pyre.applications.Script import Script as Application

# TestApp class
class TestApp(Application):
  """
  Python TestApp application.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Application.Inventory):
    """
    Python object for managing PyLithApp facilities and properties.
    """

    import pyre.inventory

    from swigtest.Scene import Scene
    scene = pyre.inventory.facility("scene", family="scene", factory=Scene)
    scene.meta['tip'] = "Scene."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="testapp"):
    """
    Constructor.
    """
    Application.__init__(self, name)
    return


  def main(self, *args, **kwds):
    """
    Run the application.
    """
    self.scene.view()

    import numpy
    x = numpy.array([1.1, 2.2, 3.3, 4.4], dtype=numpy.float32)
    self.scene.printData(x)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Application._configure(self)
    self.scene = self.inventory.scene
    return

# ======================================================================
if __name__ == '__main__':
  app = TestApp()
  app.run()


# End of file 
