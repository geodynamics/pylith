#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file tests/petsc/TestPetscApp.py
#
# @brief Test PetscApplication and PetscComponent.

import unittest

import pythia.pyre.inventory

from pylith.utils.PetscComponent import PetscComponent
from pylith.apps.PetscApplication import PetscApplication


class FooBar(PetscComponent):

    value = pythia.pyre.inventory.int("value", default=0)

    def __init__(self, name="foobar", facility="foo"):
        PetscComponent.__init__(self, name, facility)
        return

    def _cleanup(self):
        self.value = -1  # Change value to verify cleanup
        return


class TestApp(PetscApplication):
    """PetscApplication with one facility, 'foo'.
    """
    foo = pythia.pyre.inventory.facility("foo", factory=FooBar)

    def __init__(self):
        PetscApplication.__init__(self, name="testapp")
        return

    def main(self, *args, **kwds):
        self.label = "main"  # Set label to some value
        return

    def _cleanup(self):
        self.label = "_cleanup"  # Verify cleanup by changing value
        return

    def onComputeNodes(self, *args, **kwds):
        PetscApplication.onComputeNodes(self, args, kwds)
        assert(-1 == self.foo.value)
        assert("_cleanup" == self.label)
        return


class TestPetscApp(unittest.TestCase):
    """Test of PetscApplication.
    """

    def test_run(self):
        """Test running application.
        """
        TestApp().run()
        return


# ----------------------------------------------------------------------
if __name__ == '__main__':

    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestPetscApp))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
