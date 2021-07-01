# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ======================================================================
#
# @file pylith/testing/UnitTestApp.py
#
# @brief Python application for Python unit tests.

from pythia.pyre.applications.Script import Script

import unittest


def configureComponent(component):
    """Configure component and its subcomponents."""
    for subcomponent in component.components():
        configureComponent(subcomponent)
    component._configure()
    return

class TestAbstractComponent(unittest.TestCase):
    """Unit testing of abstract Pyre component.
    """
    _class = None

    def test_constructor(self):
        obj = self._class()
        self.assertTrue(obj)

    def test_configure(self):
        obj = self._class()
        self.customizeInventory(obj)
        configureComponent(obj)

    @staticmethod
    def customizeInventory(obj):
        """Customize inventory before running configure.
        """
        return


class TestComponent(TestAbstractComponent):

    _factory = None

    def test_factory(self):
        factory = self.__class__.__dict__["_factory"]
        obj = factory()
        self.assertTrue(isinstance(obj, self._class))


class UnitTestApp(Script):
    """Test application.
    """
    cov = None
    try:
        import coverage
        cov = coverage.Coverage(source=["pylith"])
    except ImportError:
        pass

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="unittestapp", petsc_options=[("malloc_dump", "true")]):
        """Constructor.
        """
        Script.__init__(self, name)
        self.petscOptions = petsc_options
        return

    def main(self):
        """Run the application.
        """
        if self.cov:
            self.cov.start()

        from pylith.utils.PetscManager import PetscManager
        petsc = PetscManager()
        petsc.options = self.petscOptions
        petsc.initialize()

        success = unittest.TextTestRunner(
            verbosity=2).run(self._suite()).wasSuccessful()

        petsc.finalize()

        if self.cov:
            self.cov.stop()
            self.cov.save()

        if not success:
            import sys
            sys.exit(1)
        return


# End of file
