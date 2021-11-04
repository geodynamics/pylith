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

from pylith.testing.FullTestApp import (FullTestCase, Check, check_data)

import meshes
import simpleshear_soln
import rigidsliding_soln

# ----------------------------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):
    
    def setUp(self):
        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": simpleshear_soln.AnalyticalSoln(),
            "mesh": self.mesh,
        }
        self.checks = [
            # #check simpleshear
            Check(
                mesh_entities=["domain", "bc_ypos"],
                vertex_fields=["displacement"],
                filename="output/simpleshear_{name}-{mesh_entity}.h5",
                exact_soln=simpleshear_soln.AnalyticalSoln(),
                defaults=defaults,
            ),
            Check(
                mesh_entities=["elastic_xpos", "elastic_xneg"],
                filename="output/simpleshear_{name}-{mesh_entity}_info.h5",
                cell_fields = ["density", "bulk_modulus", "shear_modulus"],
                exact_soln=simpleshear_soln.AnalyticalSoln(),
                defaults=defaults,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_xpos"],
                filename="output/simpleshear_{name}-{mesh_entity}_info.h5",
                vertex_fields=["initial_amplitude"],
                exact_soln=simpleshear_soln.AnalyticalSoln(),
                defaults=defaults,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_xpos"],
                filename="output/simpleshear_{name}-{mesh_entity}.h5",
                vertex_fields=["displacement"],
                exact_soln=simpleshear_soln.AnalyticalSoln(),
                defaults=defaults,
            ),
                        
            #check rigidsliding
            Check(
                mesh_entities=["domain", "bc_ypos"],
                vertex_fields=["displacement"],
                filename="output/rigidsliding_{name}-{mesh_entity}.h5",
                exact_soln=rigidsliding_soln.AnalyticalSoln(),
                defaults=defaults,
            ),
            Check(
                mesh_entities=["elastic_xpos", "elastic_xneg"],
                filename="output/rigidsliding_{name}-{mesh_entity}_info.h5",
                cell_fields = ["density", "bulk_modulus", "shear_modulus"],
                exact_soln=rigidsliding_soln.AnalyticalSoln(),
                defaults=defaults,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_xpos"],
                filename="output/rigidsliding_{name}-{mesh_entity}_info.h5",
                vertex_fields=["initial_amplitude"],
                exact_soln=rigidsliding_soln.AnalyticalSoln(),
                defaults=defaults,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_xpos"],
                filename="output/rigidsliding_{name}-{mesh_entity}.h5",
                vertex_fields=["displacement"],
                exact_soln=rigidsliding_soln.AnalyticalSoln(),
                defaults=defaults,
            ),
        ]

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args)

# ----------------------------------------------------------------------------------------------------------------------
class TestQuad(TestCase):

    def setUp(self):
        #self.name = "{problem}_quad"
        self.name = "quad"
        self.mesh = meshes.Quad
        super().setUp()

        TestCase.run_pylith(self, self.name, ["independent.cfg", "independent_quad.cfg"])
        return


# ----------------------------------------------------------------------------------------------------------------------
class TestTri(TestCase):

    def setUp(self):
        #self.name = "{problem}_tri"
        self.name = "tri"
        self.mesh = meshes.Tri
        super().setUp()

        TestCase.run_pylith(self, self.name, ["independent.cfg", "independent_tri.cfg"])
        return


# ----------------------------------------------------------------------------------------------------------------------
def test_cases():
    return [
        TestQuad
    ]

# def test_cases():
#     return [
#         TestQuad,
#         TestTri,
#     ]


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    FullTestCase.parse_args()

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
