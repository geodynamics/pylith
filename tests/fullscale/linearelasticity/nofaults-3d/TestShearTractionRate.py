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
# @file tests/fullscale/linearelasticity/nofaults-3d/TestShearTractionRate.py
#
# @brief Test suite for testing pylith with 3-D time-dependent simple shear.

import unittest

from pylith.testing.FullTestApp import check_data
from pylith.testing.FullTestApp import TestCase as FullTestCase

import meshes
from sheartraction_rate_soln import AnalyticalSoln
from sheartraction_rate_gendb import GenerateDB


# ----------------------------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):
    """Test suite for testing PyLith with 2-D time-dependent simple shear.
    """
    DIRICHLET_BOUNDARIES = ["bc_xneg", "bc_yneg", "bc_zneg"]
    NEUMANN_BOUNDARIES = ["bc_xpos", "bc_ypos"]
    OUTPUT_BOUNDARIES = ["groundsurf"]

    def setUp(self):
        """Setup for test.
        """
        FullTestCase.setUp(self)
        self.exactsoln = AnalyticalSoln()

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args, GenerateDB)

    def test_domain_solution(self):
        filename = "output/{}-domain.h5".format(self.NAME)
        vertexFields = ["displacement"]
        check_data(filename, self, self.DOMAIN, vertexFields=vertexFields)

    def test_material_info(self):
        cellFields = ["density", "bulk_modulus", "shear_modulus"]
        for material in self.MATERIALS:
            filename = "output/{}-{}_info.h5".format(self.NAME, material)
            check_data(filename, self,
                       self.MATERIALS[material], cellFields=cellFields)

    def test_material_solution(self):
        vertexFields = ["displacement", "cauchy_strain", "cauchy_stress"]
        for material in self.MATERIALS:
            filename = "output/{}-{}.h5".format(self.NAME, material)
            check_data(filename, self,
                       self.MATERIALS[material], vertexFields=vertexFields)

    def test_bcdirichlet_info(self):
        vertexFields = ["initial_amplitude",
                        "rate_start_time", "rate_amplitude"]
        for bc in self.DIRICHLET_BOUNDARIES:
            self.exactsoln.key = bc
            filename = "output/{}-{}_info.h5".format(self.NAME, bc)
            check_data(filename, self,
                       self.BOUNDARIES[bc], vertexFields=vertexFields)

    def test_bcdirichlet_solution(self):
        vertexFields = ["displacement"]
        for bc in self.DIRICHLET_BOUNDARIES:
            filename = "output/{}-{}.h5".format(self.NAME, bc)
            check_data(filename, self,
                       self.BOUNDARIES[bc], vertexFields=vertexFields)

    def test_bcneumann_info(self):
        vertexFields = ["initial_amplitude",
                        "rate_start_time", "rate_amplitude"]
        for bc in self.NEUMANN_BOUNDARIES:
            self.exactsoln.key = bc
            filename = "output/{}-{}_info.h5".format(self.NAME, bc)
            check_data(filename, self,
                       self.BOUNDARIES[bc], vertexFields=vertexFields)

    def test_bcneumann_solution(self):
        vertexFields = ["displacement"]
        for bc in self.NEUMANN_BOUNDARIES:
            filename = "output/{}-{}.h5".format(self.NAME, bc)
            check_data(filename, self,
                       self.BOUNDARIES[bc], vertexFields=vertexFields)

    def test_boundary_solution(self):
        vertexFields = ["displacement"]
        for bc in self.OUTPUT_BOUNDARIES:
            filename = "output/{}-{}.h5".format(self.NAME, bc)
            check_data(filename, self,
                       self.BOUNDARIES[bc], vertexFields=vertexFields)


# ----------------------------------------------------------------------------------------------------------------------
class TestHex(TestCase, meshes.Hex):
    NAME = "sheartraction_rate_hex"

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(
            self, self.NAME, ["sheartraction_rate.cfg", "sheartraction_rate_hex.cfg"])
        return


# ----------------------------------------------------------------------------------------------------------------------
class TestTet(TestCase, meshes.Tet):
    NAME = "sheartraction_rate_tet"

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(
            self, self.NAME, ["sheartraction_rate.cfg", "sheartraction_rate_tet.cfg"])
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
