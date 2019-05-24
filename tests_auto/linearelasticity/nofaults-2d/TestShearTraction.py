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
# @file tests_auto/linearelasticity/nofaults/TestShearTraction.py
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
    Test suite for testing PyLith with 2-D axial extension.
    """
    NAME = None  # Set in child class.

    def setUp(self):
        """
        Setup for test.
        """
        self.exactsoln = AnalyticalSoln()
        return

    def run_pylith(self, testName, args):
        run_pylith(testName, args, GenerateDB)
        return

    def test_solution_domain(self):
        filename = "output/{}-domain.h5".format(self.NAME)
        vertexFields = ["displacement"]
        cellFields = []
        check_data(filename, vertexFields, cellFields, self, self.DOMAIN)
        return

    def test_material_info(self):
        return

    def test_material_solution(self):
        return


# ----------------------------------------------------------------------------------------------------------------------
class TestQuad(TestCase, meshes.Quad):
    NAME = "sheartraction_quad"

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(self, self.NAME, ["sheartraction.cfg", "quad.cfg"])
        return


# ----------------------------------------------------------------------------------------------------------------------
class TestTri(TestCase, meshes.Tri):
    NAME = "sheartraction_tri"

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(self, self.NAME, ["sheartraction.cfg", "tri.cfg"])
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
