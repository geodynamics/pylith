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
# @file tests/fullscale/linearelasticity/nofaults-3d/TestGravityRefState.py
#
# @brief Test suite for testing pylith with 3-D gravitational body forces with initial stress and no displacement.

import unittest

from pylith.testing.FullTestApp import (FullTestCase, Check, check_data)

import meshes
import gravity_refstate_soln
import gravity_refstate_gendb

# -------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):

    def setUp(self):
        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": gravity_refstate_soln.AnalyticalSoln(),
            "mesh": self.mesh,
        }
        self.checks = [
            Check(
                mesh_entities=["domain", "groundsurf", "points"],
                vertex_fields=["displacement"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["upper_crust", "lower_crust"],
                filename="output/{name}-{mesh_entity}_info.h5",
                cell_fields = ["density", "bulk_modulus", "shear_modulus",
                              "gravitational_acceleration", "reference_strain"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["upper_crust", "lower_crust"],
                vertex_fields = ["displacement"],
                cell_fields = ["cauchy_strain", "cauchy_stress"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_xpos", "bc_yneg", "bc_ypos", "bc_zneg"],
                filename="output/{name}-{mesh_entity}_info.h5",
                cell_fields=["initial_amplitude"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_xpos", "bc_yneg", "bc_ypos", "bc_zneg"],
                vertex_fields=["displacement"],
                defaults=defaults,
            ),
        ]

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args, gravity_refstate_gendb.GenerateDB)


# -------------------------------------------------------------------------------------------------
class TestHex(TestCase):

    def setUp(self):
        self.name = "gravity_refstate_hex"
        self.mesh = meshes.Hex()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["gravity_refstate.cfg", "gravity_refstate_hex.cfg"])
        return


# -------------------------------------------------------------------------------------------------
class TestTet(TestCase):

    def setUp(self):
        self.name = "gravity_refstate_tet"
        self.mesh = meshes.Tet()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["gravity_refstate.cfg", "gravity_refstate_tet.cfg"])
        return


# -------------------------------------------------------------------------------------------------
def test_cases():
    return [
        TestHex,
        TestTet,
    ]


# -------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    FullTestCase.parse_args()

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
