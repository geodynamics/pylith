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
import shearnoslip_soln
import shearnoslip_gendb


# -------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):

    def setUp(self):
        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": shearnoslip_soln.AnalyticalSoln(),
            "mesh": self.mesh,
        }
        self.checks = [
            Check(
                mesh_entities=["domain", "boundary_ypos", "points"],
                vertex_fields=["displacement"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["mat_xneg", "mat_xmid", "mat_xposypos", "mat_xposyneg"],
                filename="output/{name}-{mesh_entity}_info.h5",
                cell_fields = ["density", "bulk_modulus", "shear_modulus"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["mat_xneg", "mat_xmid", "mat_xposypos", "mat_xposyneg"],
                vertex_fields = ["displacement"],
                cell_fields = ["cauchy_strain"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["mat_xneg", "mat_xmid", "mat_xposypos", "mat_xposyneg"],
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
                mesh_entities=["bc_xneg", "bc_xpos", "bc_yneg", "bc_ypos"],
                vertex_fields=["displacement"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["fault_xmid"],
                vertex_fields=["slip", "lagrange_multiplier_fault"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["fault_xneg"],
                vertex_fields=["slip", "lagrange_multiplier_fault"],
                defaults=defaults,
            ),
        ]

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args, shearnoslip_gendb.GenerateDB)


# -------------------------------------------------------------------------------------------------
class TestQuadGmsh(TestCase):

    def setUp(self):
        self.name = "shearnoslip_quad"
        self.mesh = meshes.QuadGmsh()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["shearnoslip.cfg", "shearnoslip_quad.cfg"])
        return


# -------------------------------------------------------------------------------------------------
class TestTriGmsh(TestCase):

    def setUp(self):
        self.name = "shearnoslip_tri"
        self.mesh = meshes.TriGmsh()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["shearnoslip.cfg", "shearnoslip_tri.cfg"])
        return


# -------------------------------------------------------------------------------------------------
def test_cases():
    return [
        TestQuadGmsh,
        TestTriGmsh,
    ]


# -------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    FullTestCase.parse_args()

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
