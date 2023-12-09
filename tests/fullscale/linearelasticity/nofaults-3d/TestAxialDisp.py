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
# @file tests/fullscale/linearelasticity/nofaults/TestAxialDisp.py
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
                mesh_entities=["domain", "groundsurf", "points"],
                vertex_fields=["displacement"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["upper_crust", "lower_crust"],
                filename="output/{name}-{mesh_entity}_info.h5",
                cell_fields = ["density", "bulk_modulus", "shear_modulus"],
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
                vertex_fields=["initial_amplitude"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_xpos", "bc_yneg", "bc_ypos", "bc_zneg"],
                vertex_fields=["displacement"],
                defaults=defaults,
            ),
        ]

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args, axialdisp_gendb.GenerateDB)


# -------------------------------------------------------------------------------------------------
class TestHex(TestCase):

    def setUp(self):
        self.name = "axialdisp_hex"
        self.mesh = meshes.Hex()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["axialdisp.cfg", "axialdisp_hex.cfg"])
        return


# -------------------------------------------------------------------------------------------------
class TestTet(TestCase):

    def setUp(self):
        self.name = "axialdisp_tet"
        self.mesh = meshes.Tet()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["axialdisp.cfg", "axialdisp_tet.cfg"])
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
