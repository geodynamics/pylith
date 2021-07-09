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
# @file tests/pytests/materials/TestRheologies.py
#
# @brief Unit testing of Python bulk rheologies object.

import unittest

from pylith.testing.UnitTestApp import TestComponent

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


if __name__ == "__main__":
    suite = unittest.TestSuite()
    for cls in [
        TestIsotropicLinearElasticity, 
        TestIsotropicLinearMaxwell, 
        TestIsotropicLinearGenMaxwell, 
        TestIsotropicPowerLaw, 
        TestIsotropicLinearIncompElasticity,
        TestIsotropicLinearPoroelasticity, 
        ]:
        suite.addTest(unittest.makeSuite(cls))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
