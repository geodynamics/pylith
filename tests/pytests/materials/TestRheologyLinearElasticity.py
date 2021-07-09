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
# @file tests/pytests/materials/TestRheologyLinearElasticity.py
#
# @brief Unit testing of Python IsotropicLinearElasticity object.

import unittest

from pylith.materials.IsotropicLinearElasticity import IsotropicLinearElasticity
from pylith.materials.AuxSubfieldsIsotropicLinearElasticity import AuxSubfieldsIsotropicLinearElasticity
from pylith.materials.DerivedSubfieldsElasticity import DerivedSubfieldsElasticity
from pylith.materials.IsotropicLinearIncompElasticity import IsotropicLinearIncompElasticity
from pylith.tests.UnitTestApp import configureComponent


class TestIsotropicLinearElasticity(unittest.TestCase):
    """Unit testing of IsotropicLinearElasticity object.
    """

    def test_constructor(self):
        rheology = IsotropicLinearElasticity()

    def test_configure(self):
        rheology = IsotropicLinearElasticity()
        configureComponent(rheology)

    def test_factory(self):
        from pylith.materials.IsotropicLinearElasticity import elasticity_rheology
        rheology = elasticity_rheology()
        self.assertTrue(isinstance(rheology, IsotropicLinearElasticity))


class TestAuxSubfieldsIsotropicLinearElasticity(unittest.TestCase):
    """Unit testing of AuxSubfieldsIsotropicLinearElasticity object.
    """

    def test_constructor(self):
        material = AuxSubfieldsIsotropicLinearElasticity()

    def test_configure(self):
        material = AuxSubfieldsIsotropicLinearElasticity()
        configureComponent(material)

    def test_factory(self):
        from pylith.materials.AuxSubfieldsIsotropicLinearElasticity import auxiliary_subfields
        materialObj = auxiliary_subfields()
        self.assertTrue(isinstance(materialObj, AuxSubfieldsIsotropicLinearElasticity))


class TestDerivedSubfieldsElasticity(unittest.TestCase):
    """Unit testing of DerivedSubfieldsElasticity object.
    """

    def test_constructor(self):
        rheology = DerivedSubfieldsElasticity()

    def test_configure(self):
        rheology = DerivedSubfieldsElasticity()
        configureComponent(rheology)

    def test_factory(self):
        from pylith.materials.DerivedSubfieldsElasticity import derived_subfields
        rheology = derived_subfields()
        self.assertTrue(isinstance(rheology, DerivedSubfieldsElasticity))


class TestIsotropicLinearIncompElasticity(unittest.TestCase):
    """Unit testing of IsotropicLinearIncompElasticity object.
    """

    def test_constructor(self):
        rheology = IsotropicLinearIncompElasticity()

    def test_configure(self):
        rheology = IsotropicLinearIncompElasticity()
        configureComponent(rheology)

    def test_factory(self):
        from pylith.materials.IsotropicLinearIncompElasticity import incompressible_elasticity_rheology
        rheology = incompressible_elasticity_rheology()
        self.assertTrue(isinstance(rheology, IsotropicLinearIncompElasticity))


if __name__ == "__main__":
    suite = unittest.TestSuite()
    for cls in [
        TestIsotropicLinearElasticity, 
        TestAuxSubfieldsIsotropicLinearElasticity,
        TestDerivedSubfieldsElasticity,
        TestIsotropicLinearIncompElasticity,
        ]:
        suite.addTest(unittest.makeSuite(cls))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
