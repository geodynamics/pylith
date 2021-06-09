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
# @file tests/fullscale/linearelasticity/faults-2d/TestTwoBlocks.py
#
# @brief Test suite for testing pylith with 2-D fault shear displacement.

import unittest

from pylith.testing.FullTestApp import check_data
from pylith.testing.FullTestApp import TestCase as FullTestCase

import meshes
from twoblocks_soln import AnalyticalSoln


# ----------------------------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):
    """Test suite for testing PyLith with fault shear displacement.
    """
    FAULTS = ["fault"]
    DIRICHLET_BOUNDARIES = ["bc_xneg", "bc_xpos"]
    OUTPUT_BOUNDARIES = ["bc_ypos"]

    def setUp(self):
        """Setup for test.
        """
        FullTestCase.setUp(self)
        self.exactsoln = AnalyticalSoln()

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args)

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
        vertexFields = ["displacement"]
        for material in self.MATERIALS:
            filename = "output/{}-{}.h5".format(self.NAME, material)
            check_data(filename, self,
                       self.MATERIALS[material], vertexFields=vertexFields)

    def test_bcdirichlet_info(self):
        vertexFields = ["initial_amplitude"]
        for bc in self.DIRICHLET_BOUNDARIES:
            filename = "output/{}-{}_info.h5".format(self.NAME, bc)
            check_data(filename, self,
                       self.BOUNDARIES[bc], vertexFields=vertexFields)

    def test_bcdirichlet_solution(self):
        vertexFields = ["displacement"]
        for bc in self.DIRICHLET_BOUNDARIES:
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
class TestQuad(TestCase, meshes.Quad):
    NAME = "twoblocks_quad"

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(
            self, self.NAME, ["twoblocks.cfg", "twoblocks_quad.cfg"])
        return


# ----------------------------------------------------------------------------------------------------------------------
class TestTri(TestCase, meshes.Tri):
    NAME = "twoblocks_tri"

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(
            self, self.NAME, ["twoblocks.cfg", "twoblocks_tri.cfg"])
        return


# ----------------------------------------------------------------------------------------------------------------------
def test_cases():
    return [
        TestQuad,
        TestTri,
    ]


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    FullTestCase.parse_args()

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
