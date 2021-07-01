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
# @file tests/fullscale/viscoelasticity/nofaults-3d/TestAxialTractionMaxwell.py
#
# @brief Test suite for testing Maxwell material with 3-D axial extension (Neumann BC).

import unittest

from pylith.testing.FullTestApp import (FullTestCase, Check, check_data)

import meshes
import axialtraction_maxwell_soln
import axialtraction_maxwell_gendb


# -------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):

    def setUp(self):
        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": axialtraction_maxwell_soln.AnalyticalSoln(),
            "mesh": self.mesh,
        }
        self.checks = [
            Check(
                mesh_entities=["domain"],
                vertex_fields=["displacement"],
                defaults=defaults,
                tolerance=2.0e-5,
            ),
            Check(
                mesh_entities=["viscomat"],
                filename="output/{name}-{mesh_entity}_info.h5",
                cell_fields = ["density", "bulk_modulus", "shear_modulus", "maxwell_time"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["viscomat"],
                vertex_fields = ["displacement", "cauchy_strain", "cauchy_stress", "viscous_strain"],
                defaults=defaults,
                tolerance=2.0e-5,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_xpos", "bc_ypos", "bc_yneg", "bc_zneg", "bc_zpos"],
                filename="output/{name}-{mesh_entity}_info.h5",
                cell_fields=["initial_amplitude"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_xpos", "bc_ypos", "bc_yneg", "bc_zneg", "bc_zpos"],
                vertex_fields=["displacement"],
                defaults=defaults,
                tolerance=2.0e-5,
            ),
        ]

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args, axialtraction_maxwell_gendb.GenerateDB)


# -------------------------------------------------------------------------------------------------
class TestHex(TestCase):

    def setUp(self):
        self.name = "axialtraction_maxwell_hex"
        self.mesh = meshes.Hex()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["axialtraction_maxwell.cfg", "axialtraction_maxwell_hex.cfg"])
        return


# -------------------------------------------------------------------------------------------------
class TestTet(TestCase):

    def setUp(self):
        self.name = "axialtraction_maxwell_tet"
        self.mesh = meshes.Tet()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["axialtraction_maxwell.cfg", "axialtraction_maxwell_tet.cfg"])
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
