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
# @file tests_auto/linearelasticity/nofaults-2d/TestShearTraction.py
#
# @brief Test suite for testing pylith with 2-D simple shear.

import unittest

from pylith.tests.FullTestApp import run_pylith, check_data

import meshes
from sheartraction_soln import AnalyticalSoln
from sheartraction_gendb import GenerateDB


# ----------------------------------------------------------------------------------------------------------------------
class TestCase(unittest.TestCase):
    """
    Test suite for testing PyLith with 2-D simple shear.
    """
    NAME = None  # Set in child class.
    DIRICHLET_BOUNDARIES = ["bc_xneg", "bc_yneg"]
    NEUMANN_BOUNDARIES = ["bc_xpos", "bc_ypos"]

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

    def test_domain_solution(self):
        filename = "output/{}-domain.h5".format(self.NAME)
        vertexFields = ["displacement"]
        check_data(filename, self, self.DOMAIN, vertexFields=vertexFields)
        return

    def test_material_info(self):
        cellFields = ["density", "bulk_modulus", "shear_modulus"]
        for material in self.MATERIALS.keys():
            filename = "output/{}-{}_info.h5".format(self.NAME, material)
            check_data(filename, self, self.MATERIALS[material], cellFields=cellFields)
        return

    def test_material_solution(self):
        vertexFields = ["displacement", "cauchy_strain", "cauchy_stress"]
        for material in self.MATERIALS.keys():
            filename = "output/{}-{}.h5".format(self.NAME, material)
            check_data(filename, self, self.MATERIALS[material], vertexFields=vertexFields)
        return

    def test_bcdirichlet_info(self):
        vertexFields = ["initial_amplitude"]
        for bc in self.DIRICHLET_BOUNDARIES:
            self.exactsoln.key = bc
            filename = "output/{}-{}_info.h5".format(self.NAME, bc)
            check_data(filename, self, self.BOUNDARIES[bc], vertexFields=vertexFields)
        return

    def test_bcdirichlet_solution(self):
        vertexFields = ["displacement"]
        for bc in self.DIRICHLET_BOUNDARIES:
            filename = "output/{}-{}.h5".format(self.NAME, bc)
            check_data(filename, self, self.BOUNDARIES[bc], vertexFields=vertexFields)
        return

    def test_bcneumann_info(self):
        vertexFields = ["initial_amplitude"]
        for bc in self.NEUMANN_BOUNDARIES:
            self.exactsoln.key = bc
            filename = "output/{}-{}_info.h5".format(self.NAME, bc)
            check_data(filename, self, self.BOUNDARIES[bc], vertexFields=vertexFields)
        return

    def test_bcneumann_solution(self):
        vertexFields = ["displacement"]
        for bc in self.NEUMANN_BOUNDARIES:
            filename = "output/{}-{}.h5".format(self.NAME, bc)
            check_data(filename, self, self.BOUNDARIES[bc], vertexFields=vertexFields)
        return


# ----------------------------------------------------------------------------------------------------------------------
class TestQuad(TestCase, meshes.Quad):
    NAME = "sheartraction_quad"

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(self, self.NAME, ["sheartraction.cfg", "sheartraction_quad.cfg"])
        return


# ----------------------------------------------------------------------------------------------------------------------
class TestTri(TestCase, meshes.Tri):
    NAME = "sheartraction_tri"

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(self, self.NAME, ["sheartraction.cfg", "sheartraction_tri.cfg"])
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
