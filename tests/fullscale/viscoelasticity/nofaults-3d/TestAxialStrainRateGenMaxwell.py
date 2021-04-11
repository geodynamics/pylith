#!/usr/bin/env nemesis
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file tests/fullscale/viscoelasticity/nofaults-3d/TestAxialStrainRateGenMaxwell.py
#
# @brief Test suite for testing generalized Maxwell material with 3-D axial strain rate (Dirichlet BC).

import unittest

from pylith.testing.FullTestApp import check_data
from pylith.testing.FullTestApp import TestCase as FullTestCase

import meshes
from axialstrainrate_genmaxwell_soln import AnalyticalSoln
from axialstrainrate_genmaxwell_gendb import GenerateDB


# ----------------------------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):
    """Test suite for testing generalized Maxwell material with 3-D axial extension (Dirichlet velocity BC).
    """
    DIRICHLET_BOUNDARIES = ["bc_xneg", "bc_xpos",
                            "bc_ypos", "bc_yneg", "bc_zneg", "bc_zpos"]

    def setUp(self):
        """Setup for test.
        """
        FullTestCase.setUp(self)
        self.exactsoln = AnalyticalSoln()
        return

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args, GenerateDB)
        return

    def test_domain_solution(self):
        filename = "output/{}-domain.h5".format(self.NAME)
        vertexFields = ["displacement"]
        check_data(filename, self, self.DOMAIN, vertexFields=vertexFields)
        return

    def test_material_info(self):
        vertexFields = ["density", "bulk_modulus",
                        "shear_modulus", "shear_modulus_ratio", "maxwell_time"]
        for material in self.MATERIALS.keys():
            filename = "output/{}-{}_info.h5".format(self.NAME, material)
            check_data(filename, self,
                       self.MATERIALS[material], vertexFields=vertexFields)
        return

    def test_material_solution(self):
        vertexFields = ["displacement", "cauchy_strain",
                        "viscous_strain", "cauchy_stress"]
        for material in self.MATERIALS.keys():
            filename = "output/{}-{}.h5".format(self.NAME, material)
            check_data(filename, self,
                       self.MATERIALS[material], vertexFields=vertexFields)
        return

    def test_bcdirichlet_info(self):
        vertexFields = ["initial_amplitude",
                        "rate_start_time", "rate_amplitude"]
        for bc in self.DIRICHLET_BOUNDARIES:
            self.exactsoln.key = bc
            filename = "output/{}-{}_info.h5".format(self.NAME, bc)
            check_data(filename, self,
                       self.BOUNDARIES[bc], vertexFields=vertexFields)
        return

    def test_bcdirichlet_solution(self):
        vertexFields = ["displacement"]
        for bc in self.DIRICHLET_BOUNDARIES:
            filename = "output/{}-{}.h5".format(self.NAME, bc)
            check_data(filename, self,
                       self.BOUNDARIES[bc], vertexFields=vertexFields)
        return


# ----------------------------------------------------------------------------------------------------------------------
class TestHex(TestCase, meshes.Hex):
    NAME = "axialstrainrate_genmaxwell_hex"

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(self, self.NAME, [
                            "axialstrainrate_genmaxwell.cfg", "axialstrainrate_genmaxwell_hex.cfg"])
        return


# ----------------------------------------------------------------------------------------------------------------------
class TestTet(TestCase, meshes.Tet):
    NAME = "axialstrainrate_genmaxwell_tet"

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(self, self.NAME, [
                            "axialstrainrate_genmaxwell.cfg", "axialstrainrate_genmaxwell_tet.cfg"])
        return


# ----------------------------------------------------------------------------------------------------------------------
def test_cases():
    return [
        TestHex,
        TestTet,
    ]


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    FullTestCase.parse_args()

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
