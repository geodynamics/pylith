#!/usr/bin/env nemesis
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#
# @file tests/fullscale/linearelasticity/cylinder-2d/TestWedgeTraction.py
#
# @brief Test suite for testing pylith with 2-D tractions on a cylinder.

import unittest

from pylith.testing.FullTestApp import (FullTestCase, Check, check_data)

import meshes
import wedgetraction_soln


# -------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):
    """Test suite for testing PyLith with normal tractions on a cylinder.
    """
    def setUp(self):
        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": wedgetraction_soln.AnalyticalSoln(),
            "mesh": self.mesh,
        }
        self.checks = [
            Check(
                mesh_entities=["domain"],
                vertex_fields=["displacement"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["elastic"],
                filename="output/{name}-{mesh_entity}_info.h5",
                cell_fields = ["density", "bulk_modulus", "shear_modulus"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["elastic"],
                vertex_fields = ["displacement", "cauchy_strain", "cauchy_stress"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_yneg", "bc_outer"],
                filename="output/{name}-{mesh_entity}_info.h5",
                vertex_fields=["initial_amplitude"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_yneg", "bc_outer"],
                vertex_fields=["displacement"],
                defaults=defaults,
            ),
        ]

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args)


# -------------------------------------------------------------------------------------------------
class TestB1Q1Quad(TestCase):

    def setUp(self):
        self.name = "wedgetraction_b1q1_quad"
        self.mesh = meshes.Quad()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["wedgetraction.cfg", "wedgetraction_b1q1_quad.cfg"])
        return


# -------------------------------------------------------------------------------------------------
class TestB2Q2Quad(TestCase):

    def setUp(self):
        self.name = "wedgetraction_b2q2_quad"
        self.mesh = meshes.Quad()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["wedgetraction.cfg", "wedgetraction_b2q2_quad.cfg"])
        return


# -------------------------------------------------------------------------------------------------
class TestB1Q1Tri(TestCase):

    def setUp(self):
        self.name = "wedgetraction_b1q1_tri"
        self.mesh = meshes.Tri()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["wedgetraction.cfg", "wedgetraction_b1q1_tri.cfg"])
        return


# -------------------------------------------------------------------------------------------------
class TestB2Q2Tri(TestCase):

    def setUp(self):
        self.name = "wedgetraction_b2q2_tri"
        self.mesh = meshes.Tri()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["wedgetraction.cfg", "wedgetraction_b2q2_tri.cfg"])
        return


# -------------------------------------------------------------------------------------------------
def test_cases():
    return [
        TestB1Q1Quad,
        TestB2Q2Quad,
        TestB1Q1Tri,
        TestB2Q2Tri,
    ]


# -------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    FullTestCase.parse_args()

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
