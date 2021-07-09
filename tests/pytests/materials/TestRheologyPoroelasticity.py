#!/usr/bin/env nemesis
#
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
