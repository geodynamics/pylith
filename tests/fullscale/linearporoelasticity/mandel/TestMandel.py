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
# @file tests/fullscale/linearporoelasticity/mandel/TestMandel.py
#
# @brief Test suite for testing pylith with Mandel's problem.

import unittest

from pylith.tests.FullTestApp import check_data
from pylith.tests.FullTestApp import TestCase as FullTestCase

import meshes
from mandel_soln import AnalyticalSoln
from mandel_gendb import GenerateDB

ratio_tolerance = 1e-0
diff_tolerance = 1e-2
# ----------------------------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):
    """
    Test suite for testing PyLith with one dimensional poroelasticity
    by means of Mandel's problem.
    """
    DIRICHLET_BOUNDARIES = ["x_neg", "x_pos", "y_neg", "y_pos"]
    #NEUMANN_BOUNDARIES = ["y_pos"]

    def setUp(self):
        """
        Setup for test.
        """
        FullTestCase.setUp(self)
        self.exactsoln = AnalyticalSoln()
        return

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args, GenerateDB)
        return

    def test_domain_solution(self):
        filename = "output/{}-domain.h5".format(self.NAME)
        vertexFields = ["displacement", "pressure", "trace_strain"]
        check_data(filename, self, self.DOMAIN, vertexFields=vertexFields, ratio_tolerance=ratio_tolerance, diff_tolerance=diff_tolerance)
        return

    def test_material_info(self):
        vertexFields = ["solid_density", "fluid_density", "fluid_viscosity", "biot_modulus", "shear_modulus", "drained_bulk_modulus", "biot_coefficient", "isotropic_permeability"]
        for material in self.MATERIALS.keys():
            filename = "output/{}-{}_info.h5".format(self.NAME, material)
            check_data(filename, self, self.MATERIALS[material], vertexFields=vertexFields, ratio_tolerance=ratio_tolerance, diff_tolerance=diff_tolerance)
        return

    def test_material_solution(self):
        vertexFields = ["displacement", "pressure", "trace_strain"]
        for material in self.MATERIALS.keys():
            filename = "output/{}-{}.h5".format(self.NAME, material)
            check_data(filename, self, self.MATERIALS[material], vertexFields=vertexFields, ratio_tolerance=ratio_tolerance, diff_tolerance=diff_tolerance)
        return

    def test_bcdirichlet_info(self):
        vertexFields = ["initial_amplitude"]
        for bc in self.DIRICHLET_BOUNDARIES:
            self.exactsoln.key = bc
            filename = "output/{}-{}_info.h5".format(self.NAME, bc)
            check_data(filename, self, self.BOUNDARIES[bc], vertexFields=vertexFields, ratio_tolerance=ratio_tolerance, diff_tolerance=diff_tolerance)
        return

    def test_bcdirichlet_solution(self):
        vertexFields = ["displacement", "pressure", "trace_strain"]
        for bc in self.DIRICHLET_BOUNDARIES:
            filename = "output/{}-{}.h5".format(self.NAME, bc)
            check_data(filename, self, self.BOUNDARIES[bc], vertexFields=vertexFields, ratio_tolerance=ratio_tolerance, diff_tolerance=diff_tolerance)
        return

# ----------------------------------------------------------------------------------------------------------------------
class TestQuad(TestCase, meshes.Quad):
    NAME = "mandel_quad"

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(self, self.NAME, ["mandel.cfg", "mandel_quad.cfg"])
        return


# ----------------------------------------------------------------------------------------------------------------------
class TestTri(TestCase, meshes.Tri):
    NAME = "mandel_tri"

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(self, self.NAME, ["mandel.cfg", "mandel_tri.cfg"])
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
