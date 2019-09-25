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
# @file tests_auto/viscoelasticity/nofaults-2d/TestAxialTractionMaxwell.py
#
# @brief Test suite for testing Maxwell material with 2-D axial extension (Neumann BC).

import unittest
import pdb

from pylith.tests.FullTestApp import check_data
from pylith.tests.FullTestApp import TestCase as FullTestCase

import meshes
from axialtraction_maxwell_soln import AnalyticalSoln


# ----------------------------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):
    """
    Test suite for testing Maxwell material with 2-D axial extension (Neumann BC).
    """
    DIRICHLET_BOUNDARIES = ["bc_xneg", "bc_ypos", "bc_yneg"]
    NEUMANN_BOUNDARIES = ["bc_xpos"]

    def setUp(self):
        """
        Setup for test.
        """
        pdb.set_trace()
        FullTestCase.setUp(self)
        self.exactsoln = AnalyticalSoln()
        return

    def run_pylith(self, testName, args):
        pdb.set_trace()
        FullTestCase.run_pylith(self, testName, args)
        return

    def test_domain_solution(self):
        pdb.set_trace()
        filename = "output/{}-domain.h5".format(self.NAME)
        vertexFields = ["displacement"]
        check_data(filename, self, self.DOMAIN, vertexFields=vertexFields)
        return

    def test_material_info(self):
        pdb.set_trace()
        cellFields = ["density", "bulk_modulus", "shear_modulus", "maxwell_time"]
        for material in self.MATERIALS.keys():
            filename = "output/{}-{}_info.h5".format(self.NAME, material)
            check_data(filename, self, self.MATERIALS[material], cellFields=cellFields)
        return

    def test_material_solution(self):
        pdb.set_trace()
        vertexFields = ["displacement", "cauchy_strain", "cauchy_stress", "viscous_strain"]
        for material in self.MATERIALS.keys():
            filename = "output/{}-{}.h5".format(self.NAME, material)
            check_data(filename, self, self.MATERIALS[material], vertexFields=vertexFields)
        return

    def test_bcdirichlet_info(self):
        pdb.set_trace()
        vertexFields = ["initial_amplitude"]
        for bc in self.DIRICHLET_BOUNDARIES:
            self.exactsoln.key = bc
            filename = "output/{}-{}_info.h5".format(self.NAME, bc)
            check_data(filename, self, self.BOUNDARIES[bc], vertexFields=vertexFields)
        return

    def test_bcdirichlet_solution(self):
        pdb.set_trace()
        vertexFields = ["displacement"]
        for bc in self.DIRICHLET_BOUNDARIES:
            filename = "output/{}-{}.h5".format(self.NAME, bc)
            check_data(filename, self, self.BOUNDARIES[bc], vertexFields=vertexFields)
        return

    def test_bcneumann_info(self):
        pdb.set_trace()
        vertexFields = ["initial_amplitude"]
        for bc in self.NEUMANN_BOUNDARIES:
            self.exactsoln.key = bc
            filename = "output/{}-{}.h5".format(self.NAME, bc)
            check_data(filename, self, self.BOUNDARIES[bc], vertexFields=vertexFields)
        return

    def test_bcneumann_solution(self):
        pdb.set_trace()
        vertexFields = ["displacement"]
        for bc in self.NEUMANN_BOUNDARIES:
            filename = "output/{}-{}.h5".format(self.NAME, bc)
            check_data(filename, self, self.BOUNDARIES[bc], vertexFields=vertexFields)
        return


# ----------------------------------------------------------------------------------------------------------------------
class TestQuad(TestCase, meshes.Quad):
    NAME = "axialtraction_maxwell_quad4"

    def setUp(self):
        pdb.set_trace()
        TestCase.setUp(self)
        TestCase.run_pylith(self, self.NAME, ["axialtraction_maxwell.cfg", "axialtraction_maxwell_quad4.cfg"])
        return


# ----------------------------------------------------------------------------------------------------------------------
class TestTri(TestCase, meshes.Tri):
    NAME = "axialtraction_maxwell_tri3"

    def setUp(self):
        pdb.set_trace()
        TestCase.setUp(self)
        TestCase.run_pylith(self, self.NAME, ["axialtraction_maxwell.cfg", "axialtraction_maxwell_tri3.cfg"])
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

    pdb.set_trace()
    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
