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
# @file tests/fullscale/linearelasticity/nofaults-2d/TestGravityIncompressible.py
#
# @brief Test suite for testing pylith with 2-D gravitational body forces for incompssible elasticity.

import unittest

from pylith.testing.FullTestApp import check_data
from pylith.testing.FullTestApp import TestCase as FullTestCase

import meshes
from gravity_incompressible_soln import AnalyticalSoln


# ----------------------------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):
    """Test suite for testing PyLith with gravitational body forces (no initial stress).
    """
    DIRICHLET_BOUNDARIES = ["bc_xneg", "bc_xpos", "bc_yneg", "bc_ypos"]

    def setUp(self):
        """Setup for test.
        """
        FullTestCase.setUp(self)
        self.exactsoln = AnalyticalSoln()
        return

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args)
        return

    def test_domain_solution(self):
        filename = "output/{}-domain.h5".format(self.NAME)
        vertexFields = ["displacement", "pressure"]
        check_data(filename, self, self.DOMAIN, vertexFields=vertexFields)
        return

    def test_material_info(self):
        cellFields = ["density", "bulk_modulus",
                      "shear_modulus", "gravitational_acceleration"]
        for material in self.MATERIALS.keys():
            filename = "output/{}-{}_info.h5".format(self.NAME, material)
            check_data(filename, self,
                       self.MATERIALS[material], cellFields=cellFields)
        return

    def test_material_solution(self):
        vertexFields = ["displacement", "pressure",
                        "cauchy_strain", "cauchy_stress"]
        for material in self.MATERIALS.keys():
            filename = "output/{}-{}.h5".format(self.NAME, material)
            check_data(filename, self,
                       self.MATERIALS[material], vertexFields=vertexFields)
        return

    def test_bcdirichlet_info(self):
        cellFields = ["initial_amplitude"]
        for bc in self.DIRICHLET_BOUNDARIES:
            self.exactsoln.key = bc
            filename = "output/{}-{}_info.h5".format(self.NAME, bc)
            check_data(filename, self,
                       self.BOUNDARIES[bc], cellFields=cellFields)
        return

    def test_bcdirichlet_solution(self):
        vertexFields = ["displacement", "pressure"]
        for bc in self.DIRICHLET_BOUNDARIES:
            filename = "output/{}-{}.h5".format(self.NAME, bc)
            check_data(filename, self,
                       self.BOUNDARIES[bc], vertexFields=vertexFields)
        return


# ----------------------------------------------------------------------------------------------------------------------
class TestQuad(TestCase, meshes.Quad):
    NAME = "gravity_incompressible_quad"

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(self, self.NAME, [
                            "gravity_incompressible.cfg", "gravity_incompressible_quad.cfg"])
        return


# ----------------------------------------------------------------------------------------------------------------------
class TestTri(TestCase, meshes.Tri):
    NAME = "gravity_incompressible_tri"

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(self, self.NAME, [
                            "gravity_incompressible.cfg", "gravity_incompressible_tri.cfg"])
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
