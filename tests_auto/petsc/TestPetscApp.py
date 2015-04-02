#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file tests/petsc/TestPetscApp.py
##
## @brief Test PetscApplication and PetscComponent.

import unittest

# Component with value.
from pylith.utils.PetscComponent import PetscComponent
class FooBar(PetscComponent):
  class Inventory(PetscComponent.Inventory):
    import pyre.inventory
    value = pyre.inventory.int("value", default=0)
  def __init__(self, name="foobar", facility="foo"):
    PetscComponent.__init__(self, name, facility)
    return
  def _configure(self):
    PetscComponent._configure(self)
    self.value = self.inventory.value
    return
  def _cleanup(self):
    self.value = -1 # Change value to verify cleanup
    return


# Local version of PetscApplicatiion
from pylith.apps.PetscApplication import PetscApplication
class TestApp(PetscApplication):
  """
  PetscApplication with one facility, 'foo'.
  """
  class Inventory(PetscApplication.Inventory):
    import pyre.inventory
    foo = pyre.inventory.facility("foo", factory=FooBar)
  def __init__(self):
    PetscApplication.__init__(self, name="testapp")
    return
  def main(self, *args, **kwds):
    self.label = "main" # Set label to some value
    return
  def _configure(self):
    PetscApplication._configure(self)
    self.foo = self.inventory.foo
    return
  def _cleanup(self):
    self.label = "_cleanup" # Verify cleanup by changing value
    return
  def onComputeNodes(self, *args, **kwds):
    PetscApplication.onComputeNodes(self, args, kwds)
    assert(-1 == self.foo.value)
    assert("_cleanup" == self.label)
    return


class TestPetscApp(unittest.TestCase):
  """
  Test of PetscApplication.
  """

  def test_run(self):
    """
    Test running application.
    """
    app = TestApp()
    app.run()

    return


# ----------------------------------------------------------------------
if __name__ == '__main__':
  import unittest
  from TestPetscApp import TestPetscApp as Tester

  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(Tester))
  unittest.TextTestRunner(verbosity=2).run(suite)


# End of file 
