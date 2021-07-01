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
# @file tests/fullscale/linearelasticity/nofaults-2d/TestAxialDisp.py
#
# @brief Test suite for testing pylith with 2-D axial extension.

import unittest

from pylith.testing.FullTestApp import (FullTestCase, Check, check_data)

import meshes
import axialdisp_soln
import axialdisp_gendb


# -------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):

    def setUp(self):
        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": axialdisp_soln.AnalyticalSoln(),
            "mesh": self.mesh,
        }
        self.checks = [
            Check(
                mesh_entities=["domain", "bc_ypos", "points"],
                vertex_fields=["displacement"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["elastic_xpos", "elastic_xneg"],
                filename="output/{name}-{mesh_entity}_info.h5",
                cell_fields = ["density", "bulk_modulus", "shear_modulus"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["elastic_xpos", "elastic_xneg"],
                vertex_fields = ["displacement", "cauchy_strain", "cauchy_stress"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_xpos", "bc_yneg"],
                filename="output/{name}-{mesh_entity}_info.h5",
                vertex_fields=["initial_amplitude"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_xpos", "bc_yneg"],
                vertex_fields=["displacement"],
                defaults=defaults,
            ),
        ]

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args, axialdisp_gendb.GenerateDB)


# -------------------------------------------------------------------------------------------------
class TestQuad(TestCase):

    def setUp(self):
        self.name = "axialdisp_quad"
        self.mesh = meshes.Quad()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["axialdisp.cfg", "axialdisp_quad.cfg"])
        return


# -------------------------------------------------------------------------------------------------
class TestTri(TestCase):

    def setUp(self):
        self.name = "axialdisp_tri"
        self.mesh = meshes.Tri()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["axialdisp.cfg", "axialdisp_tri.cfg"])
        return


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
