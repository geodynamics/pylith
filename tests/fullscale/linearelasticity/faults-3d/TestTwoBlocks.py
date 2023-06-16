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
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

import unittest

from pylith.testing.FullTestApp import (FullTestCase, Check)

import meshes
import twoblocks_soln


# -------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):

    def setUp(self):
        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": twoblocks_soln.AnalyticalSoln(),
            "mesh": self.mesh,
        }
        self.checks = [
            Check(
                mesh_entities=["domain", "points"],
                vertex_fields=["displacement"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["mat_elastic"],
                filename="output/{name}-{mesh_entity}_info.h5",
                cell_fields = ["density", "bulk_modulus", "shear_modulus"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["mat_elastic"],
                vertex_fields = ["displacement"],
                cell_fields = ["cauchy_strain"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["mat_elastic"],
                cell_fields = ["cauchy_stress"],
                scale = 1.0e+6,
                defaults=defaults,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_xpos"],
                filename="output/{name}-{mesh_entity}_info.h5",
                vertex_fields=["initial_amplitude"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_xpos"],
                vertex_fields=["displacement"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["fault"],
                vertex_fields=["slip"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["fault"],
                vertex_fields=["lagrange_multiplier_fault"],
                scale = 1.0e+6,
                defaults=defaults,
            ),
        ]

    def run_pylith(self, testName, args, nprocs=1):
        FullTestCase.run_pylith(self, testName, args, nprocs=nprocs)


# -------------------------------------------------------------------------------------------------
class TestHexGmsh(TestCase):

    def setUp(self):
        self.name = "twoblocks_hex"
        self.mesh = meshes.HexGmsh()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["twoblocks.cfg", "twoblocks_hex.cfg"], nprocs=2)
        return


# -------------------------------------------------------------------------------------------------
class TestTetGmsh(TestCase):

    def setUp(self):
        self.name = "twoblocks_tet"
        self.mesh = meshes.TetGmsh()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["twoblocks.cfg", "twoblocks_tet.cfg"], nprocs=3)
        return


# -------------------------------------------------------------------------------------------------
def test_cases():
    return [
        TestHexGmsh,
        TestTetGmsh,
    ]


# -------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    FullTestCase.parse_args()

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
