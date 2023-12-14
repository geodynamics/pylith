#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file tests/pytests/materials/TestRheologyPoroelasticity.py
#
# @brief Unit testing of Python poroelasticity rheology object.

import unittest

from pylith.materials.IsotropicLinearPoroelasticity import IsotropicLinearPoroelasticity
from pylith.tests.UnitTestApp import configureComponent


class TestIsotropicLinearPoroelasticity(unittest.TestCase):
    """Unit testing of IsotropicLinearPoroelasticity object.
    """

    def test_constructor(self):
        rheology = IsotropicLinearPoroelasticity()

    def test_configure(self):
        rheology = IsotropicLinearPoroelasticity()
        configureComponent(rheology)

    def test_factory(self):
        from pylith.materials.IsotropicLinearPoroelasticity import poroelasticity_rheology
        rheology = poroelasticity_rheology()
        self.assertTrue(isinstance(rheology, IsotropicLinearPoroelasticity))


if __name__ == "__main__":
    suite = unittest.TestSuite()
    for cls in [
        TestIsotropicLinearPoroelasticity,
        ]:
        suite.addTest(unittest.makeSuite(cls))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
