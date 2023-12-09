#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file tests/fullscale/poroelasticity/terzaghi_compaction/TestTerzaghiCompaction.py
#
# @brief Test suite for testing pylith with Terzaghi's problem, involving compaction (change of porosity).

import unittest

from pylith.testing.FullTestApp import (FullTestCase, Check, check_data)

import meshes
import terzaghi_compaction_soln
import terzaghi_compaction_gendb

# We do not include trace_strain in the solution fields, because of the
# poor convergence of the series solution.
SOLUTION_FIELDS = ["displacement", "pressure"]
SOLUTION_TOLERANCE = 0.2

# -------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):

    def setUp(self):
        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": terzaghi_compaction_soln.AnalyticalSoln(),
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
                cell_fields=[
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
                mesh_entities=["x_neg", "x_pos", "y_pos_dir", "y_neg", "y_pos_neu"],
                filename="output/{name}-{mesh_entity}_info.h5",
                vertex_fields=["initial_amplitude"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["x_neg", "x_pos", "y_pos_dir", "y_neg", "y_pos_neu"],
                vertex_fields=SOLUTION_FIELDS,
                defaults=defaults,
                tolerance=SOLUTION_TOLERANCE,
            ),
        ]

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args, terzaghi_compaction_gendb.GenerateDB)


# -------------------------------------------------------------------------------------------------
class TestQuad(TestCase):

    def setUp(self):
        self.name = "terzaghi_compaction_quad"
        self.mesh = meshes.Quad()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["terzaghi_compaction.cfg", "terzaghi_compaction_quad.cfg"])


# -------------------------------------------------------------------------------------------------
class TestTri(TestCase):

    def setUp(self):
        self.name = "terzaghi_compaction_tri"
        self.mesh = meshes.Tri()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["terzaghi_compaction.cfg", "terzaghi_compaction_tri.cfg"])


# -------------------------------------------------------------------------------------------------
def test_cases():
    return [
        TestQuad,
        TestTri,
    ]


# -------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    FullTestCase.parse_args()

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
