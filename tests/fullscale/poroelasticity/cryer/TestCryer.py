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
# @file tests/fullscale/poroelasticity/cryer/TestCryer.py
#
# @brief Test suite for testing pylith with Cryer's problem.

import unittest

from pylith.testing.FullTestApp import (FullTestCase, Check, check_data)

import meshes
import cryer_soln

# We do not include trace_strain in the test of the solution fields, because of the
# poor convergence of the series solution.
SOLUTION_FIELDS = ["displacement", "pressure"]
SOLUTION_TOLERANCE = 0.5

# -------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):

    def setUp(self):
        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": cryer_soln.AnalyticalSoln(),
            "mesh": self.mesh,
        }
        self.checks = [
            Check(
                mesh_entities=["domain"],
                vertex_fields=SOLUTION_FIELDS,
                defaults=defaults,
                tolerance=SOLUTION_TOLERANCE,
            ),
            Check(
                mesh_entities=["poroelastic"],
                filename="output/{name}-{mesh_entity}_info.h5",
                cell_fields = [
                    "biot_coefficient",
                    "biot_modulus",
                    "drained_bulk_modulus",
                    "fluid_density",
                    "fluid_viscosity",
                    "isotropic_permeability",
                    "porosity",
                    "shear_modulus",
                    "solid_density",
                ],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["poroelastic"],
                vertex_fields = SOLUTION_FIELDS,
                defaults=defaults,
                tolerance=SOLUTION_TOLERANCE,
            ),
            Check(
                mesh_entities=["x_neg", "y_neg", "z_neg", "surface_pressure"],
                filename="output/{name}-{mesh_entity}_info.h5",
                vertex_fields=["initial_amplitude"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["x_neg", "y_neg", "z_neg", "surface_pressure"],
                vertex_fields=SOLUTION_FIELDS,
                defaults=defaults,
                tolerance=SOLUTION_TOLERANCE,
            ),
        ]

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args)


# -------------------------------------------------------------------------------------------------
class TestHex(TestCase):

    def setUp(self):
        self.name = "cryer_hex"
        self.mesh = meshes.Hex()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["cryer.cfg", "cryer_hex.cfg"])
        return


# -------------------------------------------------------------------------------------------------
class TestTet(TestCase):

    def setUp(self):
        self.name = "cryer_tet"
        self.mesh = meshes.Tet()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["cryer.cfg", "cryer_tet.cfg"])
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
