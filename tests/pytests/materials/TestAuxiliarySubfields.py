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
# @file tests/pytests/materials/TestAuxiliarySubfields.py
#
# @brief Unit testing of Python TestAuxiliarySubfields object.

import unittest

from pylith.testing.TestCases import TestComponent, make_suite
import pylith.materials.AuxSubfieldsElasticity
import pylith.materials.AuxSubfieldsIsotropicLinearElasticity
import pylith.materials.AuxSubfieldsIsotropicLinearMaxwell
import pylith.materials.AuxSubfieldsIsotropicLinearGenMaxwell
import pylith.materials.AuxSubfieldsIsotropicPowerLaw
import pylith.materials.AuxSubfieldsPoroelasticity
import pylith.materials.AuxSubfieldsIsotropicLinearPoroelasticity


class TestAuxSubfieldsElasticity(TestComponent):
    """Unit testing of AuxSubfieldsElasticity object.
    """
    _class = pylith.materials.AuxSubfieldsElasticity.AuxSubfieldsElasticity
    _factory = pylith.materials.AuxSubfieldsElasticity.auxiliary_subfields


class TestAuxSubfieldsIsotropicLinearElasticity(TestComponent):
    """Unit testing of AuxSubfieldsIsotropicLinearElasticity object.
    """
    _class = pylith.materials.AuxSubfieldsIsotropicLinearElasticity.AuxSubfieldsIsotropicLinearElasticity
    _factory = pylith.materials.AuxSubfieldsIsotropicLinearElasticity.auxiliary_subfields


class TestAuxSubfieldsIsotropicLinearMaxwell(TestComponent):
    """Unit testing of AuxSubfieldsIsotropicLinearMaxwell object.
    """
    _class = pylith.materials.AuxSubfieldsIsotropicLinearMaxwell.AuxSubfieldsIsotropicLinearMaxwell
    _factory = pylith.materials.AuxSubfieldsIsotropicLinearMaxwell.auxiliary_subfields


class TestAuxSubfieldsIsotropicLinearGenMaxwell(TestComponent):
    """Unit testing of AuxSubfieldsIsotropicLinearGenMaxwell object.
    """
    _class = pylith.materials.AuxSubfieldsIsotropicLinearGenMaxwell.AuxSubfieldsIsotropicLinearGenMaxwell
    _factory = pylith.materials.AuxSubfieldsIsotropicLinearGenMaxwell.auxiliary_subfields


class TestAuxSubfieldsIsotropicPowerLaw(TestComponent):
    """Unit testing of AuxSubfieldsIsotropicPowerLaw object.
    """
    _class = pylith.materials.AuxSubfieldsIsotropicPowerLaw.AuxSubfieldsIsotropicPowerLaw
    _factory = pylith.materials.AuxSubfieldsIsotropicPowerLaw.auxiliary_subfields


class TestAuxSubfieldsPoroelasticity(TestComponent):
    """Unit testing of AuxSubfieldsPoroelasticity object.
    """
    _class = pylith.materials.AuxSubfieldsPoroelasticity.AuxSubfieldsPoroelasticity
    _factory = pylith.materials.AuxSubfieldsPoroelasticity.auxiliary_subfields


class TestAuxSubfieldsIsotropicLinearPoroelasticity(TestComponent):
    """Unit testing of AuxSubfieldsIsotropicLinearPoroelasticity object.
    """
    _class = pylith.materials.AuxSubfieldsIsotropicLinearPoroelasticity.AuxSubfieldsIsotropicLinearPoroelasticity
    _factory = pylith.materials.AuxSubfieldsIsotropicLinearPoroelasticity.auxiliary_subfields


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [
        TestAuxSubfieldsElasticity,
        TestAuxSubfieldsIsotropicLinearElasticity,
        TestAuxSubfieldsIsotropicLinearMaxwell,
        TestAuxSubfieldsIsotropicLinearGenMaxwell,
        TestAuxSubfieldsIsotropicPowerLaw,
        TestAuxSubfieldsPoroelasticity,
        TestAuxSubfieldsIsotropicLinearPoroelasticity,
        ]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
