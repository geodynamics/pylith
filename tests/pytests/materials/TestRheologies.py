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

from pylith.testing.TestCases import TestComponent, make_suite

import pylith.materials.IsotropicLinearElasticity
import pylith.materials.IsotropicLinearMaxwell
import pylith.materials.IsotropicLinearGenMaxwell
import pylith.materials.IsotropicPowerLaw

import pylith.materials.IsotropicLinearIncompElasticity

import pylith.materials.IsotropicLinearPoroelasticity


class TestIsotropicLinearElasticity(TestComponent):
    """Unit testing of IsotropicLinearElasticity object.
    """
    _class = pylith.materials.IsotropicLinearElasticity.IsotropicLinearElasticity
    _factory = pylith.materials.IsotropicLinearElasticity.elasticity_rheology


class TestIsotropicLinearMaxwell(TestComponent):
    """Unit testing of IsotropicLinearMaxwell object.
    """
    _class = pylith.materials.IsotropicLinearMaxwell.IsotropicLinearMaxwell
    _factory = pylith.materials.IsotropicLinearMaxwell.elasticity_rheology


class TestIsotropicLinearGenMaxwell(TestComponent):
    """Unit testing of IsotropicLinearGenMaxwell object.
    """
    _class = pylith.materials.IsotropicLinearGenMaxwell.IsotropicLinearGenMaxwell
    _factory = pylith.materials.IsotropicLinearGenMaxwell.elasticity_rheology


class TestIsotropicPowerLaw(TestComponent):
    """Unit testing of IsotropicPowerLaw object.
    """
    _class = pylith.materials.IsotropicPowerLaw.IsotropicPowerLaw
    _factory = pylith.materials.IsotropicPowerLaw.elasticity_rheology


class TestIsotropicLinearIncompElasticity(TestComponent):
    """Unit testing of IsotropicLinearIncompElasticity object.
    """
    _class = pylith.materials.IsotropicLinearIncompElasticity.IsotropicLinearIncompElasticity
    _factory = pylith.materials.IsotropicLinearIncompElasticity.incompressible_elasticity_rheology


class TestIsotropicLinearPoroelasticity(TestComponent):
    """Unit testing of IsotropicLinearPoroelasticity object.
    """
    _class = pylith.materials.IsotropicLinearPoroelasticity.IsotropicLinearPoroelasticity
    _factory = pylith.materials.IsotropicLinearPoroelasticity.poroelasticity_rheology


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [
        TestIsotropicLinearElasticity, 
        TestIsotropicLinearMaxwell, 
        TestIsotropicLinearGenMaxwell, 
        TestIsotropicPowerLaw, 
        TestIsotropicLinearIncompElasticity,
        TestIsotropicLinearPoroelasticity, 
    ]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
