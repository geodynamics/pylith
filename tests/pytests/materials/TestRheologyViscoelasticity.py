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
# @file tests/pytests/materials/TestRheologyViscoelasticity.py
#
# @brief Unit testing of Python viscoelasticity rheology objects.

import unittest

from pylith.materials.IsotropicLinearMaxwell import IsotropicLinearMaxwell
from pylith.materials.AuxSubfieldsIsotropicLinearMaxwell import AuxSubfieldsIsotropicLinearMaxwell
from pylith.materials.IsotropicLinearGenMaxwell import IsotropicLinearGenMaxwell
from pylith.materials.AuxSubfieldsIsotropicLinearGenMaxwell import AuxSubfieldsIsotropicLinearGenMaxwell
from pylith.materials.IsotropicPowerLaw import IsotropicPowerLaw
from pylith.materials.AuxSubfieldsIsotropicPowerLaw import AuxSubfieldsIsotropicPowerLaw
from pylith.tests.UnitTestApp import configureComponent


class TestIsotropicLinearMaxwell(unittest.TestCase):
    """Unit testing of IsotropicLinearMaxwell object.
    """

    def test_constructor(self):
        rheology = IsotropicLinearMaxwell()

    def test_configure(self):
        rheology = IsotropicLinearMaxwell()
        configureComponent(rheology)

    def test_factory(self):
        from pylith.materials.IsotropicLinearMaxwell import elasticity_rheology
        rheology = elasticity_rheology()
        self.assertTrue(isinstance(rheology, IsotropicLinearMaxwell))


class TestAuxSubfieldsIsotropicLinearMaxwell(unittest.TestCase):
    """Unit testing of AuxSubfieldsIsotropicLinearMaxwell object.
    """

    def test_constructor(self):
        material = AuxSubfieldsIsotropicLinearMaxwell()

    def test_configure(self):
        material = AuxSubfieldsIsotropicLinearMaxwell()
        configureComponent(material)

    def test_factory(self):
        from pylith.materials.AuxSubfieldsIsotropicLinearMaxwell import auxiliary_subfields
        materialObj = auxiliary_subfields()
        self.assertTrue(isinstance(materialObj, AuxSubfieldsIsotropicLinearMaxwell))


class TestIsotropicLinearGenMaxwell(unittest.TestCase):
    """Unit testing of IsotropicLinearGenMaxwell object.
    """

    def test_constructor(self):
        rheology = IsotropicLinearGenMaxwell()

    def test_configure(self):
        rheology = IsotropicLinearGenMaxwell()
        configureComponent(rheology)

    def test_factory(self):
        from pylith.materials.IsotropicLinearGenMaxwell import elasticity_rheology
        rheology = elasticity_rheology()
        self.assertTrue(isinstance(rheology, IsotropicLinearGenMaxwell))


class TestAuxSubfieldsIsotropicLinearGenMaxwell(unittest.TestCase):
    """Unit testing of AuxSubfieldsIsotropicLinearGenMaxwell object.
    """

    def test_constructor(self):
        material = AuxSubfieldsIsotropicLinearGenMaxwell()

    def test_configure(self):
        material = AuxSubfieldsIsotropicLinearGenMaxwell()
        configureComponent(material)

    def test_factory(self):
        from pylith.materials.AuxSubfieldsIsotropicLinearGenMaxwell import auxiliary_subfields
        materialObj = auxiliary_subfields()
        self.assertTrue(isinstance(materialObj, AuxSubfieldsIsotropicLinearGenMaxwell))


class TestIsotropicPowerLaw(unittest.TestCase):
    """Unit testing of IsotropicPowerLaw object.
    """

    def test_constructor(self):
        rheology = IsotropicPowerLaw()

    def test_configure(self):
        rheology = IsotropicPowerLaw()
        configureComponent(rheology)

    def test_factory(self):
        from pylith.materials.IsotropicPowerLaw import elasticity_rheology
        rheology = elasticity_rheology()
        self.assertTrue(isinstance(rheology, IsotropicPowerLaw))


class TestAuxSubfieldsIsotropicPowerLaw(unittest.TestCase):
    """Unit testing of AuxSubfieldsIsotropicPowerLaw object.
    """

    def test_constructor(self):
        material = AuxSubfieldsIsotropicPowerLaw()

    def test_configure(self):
        material = AuxSubfieldsIsotropicPowerLaw()
        configureComponent(material)

    def test_factory(self):
        from pylith.materials.AuxSubfieldsIsotropicPowerLaw import auxiliary_subfields
        materialObj = auxiliary_subfields()
        self.assertTrue(isinstance(materialObj, AuxSubfieldsIsotropicPowerLaw))


if __name__ == "__main__":
    suite = unittest.TestSuite()
    for cls in [
        TestIsotropicLinearMaxwell,
        TestAuxSubfieldsIsotropicLinearMaxwell,
        TestIsotropicLinearGenMaxwell,
        TestAuxSubfieldsIsotropicLinearGenMaxwell,
        TestIsotropicPowerLaw,
        TestAuxSubfieldsIsotropicPowerLaw,
        ]:
        suite.addTest(unittest.makeSuite(cls))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
