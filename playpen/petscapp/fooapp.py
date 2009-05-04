#!/usr/bin/env nemesis

from pyre.applications.Script import Script as Application
from pyre.components.Component import Component

# ======================================================================
# PetscApplication
class PetscApplication(Application):
  """
  PETSc application.
  """

  def __init__(self, name="petscapp"):
    """
    Constructor.
    """
    Application.__init__(self)
    return


  def main(self, *args, **kwds):
    """
    Run simulation in parallel on compute nodes. Would be
    onComputeNodes() in a real Pyre mpi application.
    """
    print "Call PetscInitialize()."
    self.fake_main(*args, **kwds)
    self.cleanup()
    print "Call PetscFinalize()."
    return


  def cleanup(self):
    """
    Deallocate data structures.
    """
    for component in self.components():
      if isinstance(component, PetscComponent):
        component.cleanup()
    self._cleanup()
    return


  def _cleanup(self):
    """
    Deallocate locally managed data structures.
    """
    return
    

# ======================================================================
# PetscComponent
class PetscComponent(Component):
  """
  PETSc component.
  """

  def __init__(self, name="petsccomponent", facility="petsccomponent"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility)
    return


  def cleanup(self):
    """
    Deallocate data structures.
    """
    for component in self.components():
      if isinstance(component, PetscComponent):
        component.cleanup()
    self._cleanup()
    return


  def _cleanup(self):
    """
    Deallocate locally managed data structures.
    """
    return
    

# ======================================================================
# Bar
class Bar(PetscComponent):
  """
  Bar Pyre component.
  """
  
  # Inventory
  import pyre.inventory
  value = pyre.inventory.int("value", default=0)
  value.meta['tip'] = "An integer value."


  def __init__(self, name="bar"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="bar")
    return


  def show(self):
    """
    Print value.
    """
    print "  Bar value is %d" % self.value
    return


  def _cleanup(self):
    """
    Deallocate locally managed data structures.
    """
    print "Deallocate Bar data structures in _cleanup()."
    return


# ======================================================================
# FooApp class
class FooApp(PetscApplication):
  """
  FooApp application.
  """
  
  # Inventory
  import pyre.inventory
  foo = pyre.inventory.facility("foo", factory=Bar)
  foo.meta['tip'] = "Facility foo."


  def __init__(self, name="fooapp"):
    """
    Constructor.
    """
    Application.__init__(self, name)
    self.foo2 = Bar()
    self.foo2.value = 2
    return


  def fake_main(self, *args, **kwds):
    """
    Run the application. Would be called the usual main() in a real
    application.
    """
    print "Doing some stuff in main."
    print "Foo:"
    self.foo.show()
    print "Foo 2:"
    self.foo2.show()

    return
  

  def _cleanup(self):
    """
    Deallocate locally managed data structures.
    """
    print "Deallocate application data structures in _cleanup()."
    self.foo2.cleanup()
    return
    

# ======================================================================
if __name__ == "__main__":

  app = FooApp()
  app.run()


# End of file
