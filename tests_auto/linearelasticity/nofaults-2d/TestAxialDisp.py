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
# @file tests_auto/linearelasticity/nofaults/TestAxialDisp.py
##
# @brief Test suite for testing pylith with 2-D axial extension.

import unittest

from pylith.tests.FullTestApp import run_pylith, check_data

import meshes
from axialdisp_soln import AnalyticalSoln
from axialdisp_gendb import GenerateDB


# ----------------------------------------------------------------------------------------------------------------------
class TestCase(unittest.TestCase):
    """
    Test suite for testing PyLith with 2-D axial extension.
    """
    NAME = None  # Set in child class.

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
        run_pylith(testName, args, GenerateDB)
        return

    def test_solution_domain(self):
        filename = "output/{}-domain.h5".format(self.NAME)
        vertexFields = ["displacement"]
        cellFields = []
        check_data(filename, vertexFields, cellFields, self, self.DOMAIN, self.verbosity)
        return

    def test_material_info(self):
        vertexFields = []
        cellFields = ["density", "bulk_modulus", "shear_modulus"]
        for material in self.MATERIALS.keys():
            filename = "output/{}-{}_info.h5".format(self.NAME, material)
            check_data(filename, vertexFields, cellFields, self, self.MATERIALS[material], self.verbosity)
        return

    def test_material_solution(self):
        vertexFields = ["displacement"]
        cellFields = []
        for material in self.MATERIALS.keys():
            filename = "output/{}-{}.h5".format(self.NAME, material)
            check_data(filename, vertexFields, cellFields, self, self.MATERIALS[material], self.verbosity)
        return


# ----------------------------------------------------------------------------------------------------------------------
class TestQuad(TestCase, meshes.Quad):
    NAME = "axialdisp_quad"

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(self, self.NAME, ["axialdisp.cfg", "quad.cfg"])
        return


# ----------------------------------------------------------------------------------------------------------------------
class TestTri(TestCase, meshes.Tri):
    NAME = "axialdisp_tri"

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(self, self.NAME, ["axialdisp.cfg", "tri.cfg"])
        return


# ----------------------------------------------------------------------------------------------------------------------
def test_cases():
    return [
        TestQuad,
        TestTri,
    ]


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
