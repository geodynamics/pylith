#!/usr/bin/env python
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
# @file tests_auto/linearelasticity/nofaults-2d/TestGravity.py
#
# @brief Test suite for testing pylith with 2-D gravitational body forces (no initial stress).

import unittest

from pylith.tests.FullTestApp import run_pylith, check_data

import meshes
from gravity_soln import AnalyticalSoln


# ----------------------------------------------------------------------------------------------------------------------
class TestCase(unittest.TestCase):
    """
    Test suite for testing PyLith with gravitational body forces (no initial stress).
    """
    NAME = None  # Set in child class.
    DIRICHLET_BOUNDARIES = ["bc_xneg", "bc_xpos", "bc_yneg"]

    def setUp(self):
        """
        Setup for test.
        """
        self.exactsoln = AnalyticalSoln()
        self.verbosity = 0
        return

    def run_pylith(self, testName, args):
        if self.verbosity > 0:
            print("Running Pylith with args '{}' ...".format(" ".join(args)))
        #run_pylith(testName, args)
        return

    def test_solution_domain(self):
        filename = "output/{}-domain.h5".format(self.NAME)
        vertexFields = ["displacement"]
        check_data(filename, self, self.DOMAIN, vertexFields=vertexFields)
        return

    def test_material_info(self):
        cellFields = ["density", "bulk_modulus", "shear_modulus", "gravitational_acceleration"]
        for material in self.MATERIALS.keys():
            filename = "output/{}-{}_info.h5".format(self.NAME, material)
            check_data(filename, self, self.MATERIALS[material], cellFields=cellFields)
        return

    @unittest.expectedFailure
    def test_material_solution(self):
        vertexFields = ["displacement"]
        for material in self.MATERIALS.keys():
            filename = "output/{}-{}.h5".format(self.NAME, material)
            check_data(filename, self, self.MATERIALS[material], vertexFields=vertexFields)
        return

    def test_bcdirichlet_info(self):
        cellFields = ["initial_amplitude"]
        for bc in self.DIRICHLET_BOUNDARIES:
            filename = "output/{}-{}_info.h5".format(self.NAME, bc)
            check_data(filename, self, self.BOUNDARIES[bc], cellFields=cellFields)
        return

    @unittest.expectedFailure
    def test_bcdirichlet_solution(self):
        vertexFields = ["displacement"]
        for bc in self.DIRICHLET_BOUNDARIES:
            filename = "output/{}-{}.h5".format(self.NAME, bc)
            check_data(filename, self, self.BOUNDARIES[bc], vertexFields=vertexFields)
        return


# ----------------------------------------------------------------------------------------------------------------------
class TestQuad(TestCase, meshes.Quad):
    NAME = "gravity_quad"

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(self, self.NAME, ["gravity.cfg", "gravity_quad.cfg"])
        return


# ----------------------------------------------------------------------------------------------------------------------
class TestTri(TestCase, meshes.Tri):
    NAME = "gravity_tri"

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(self, self.NAME, ["gravity.cfg", "gravity_tri.cfg"])
        return


# ----------------------------------------------------------------------------------------------------------------------
def test_cases():
    return [
        TestQuad,
        # TestTri,
    ]


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
