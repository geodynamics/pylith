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
# @file tests/pytests/materials/TestAuxiliarySubfields.py
#
# @brief Unit testing of Python TestAuxiliarySubfields object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
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


if __name__ == "__main__":
    suite = unittest.TestSuite()
    for cls in [
        TestAuxSubfieldsElasticity,
        TestAuxSubfieldsIsotropicLinearElasticity,
        TestAuxSubfieldsIsotropicLinearMaxwell,
        TestAuxSubfieldsIsotropicLinearGenMaxwell,
        TestAuxSubfieldsIsotropicPowerLaw,
        TestAuxSubfieldsPoroelasticity,
        TestAuxSubfieldsIsotropicLinearPoroelasticity,
        ]:
        suite.addTest(unittest.makeSuite(cls))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
