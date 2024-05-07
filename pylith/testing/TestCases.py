# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================


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


def make_suite(test_classes, loader=unittest.defaultTestLoader):
    suite = unittest.TestSuite()
    for cls in test_classes:
        suite.addTests(loader.loadTestsFromTestCase(cls))
    return suite


# End of file
