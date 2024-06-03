#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file tests/fullscale/poroelasticity/mandel/TestMandel.py
#
# @brief Test suite for testing pylith with Mandel's problem.
#
# We do not include trace_strain in the check of the solution fields, because of the
# poor convergence of the series solution.

import unittest

from pylith.testing.FullTestApp import (FullTestCase, Check, check_data)

import meshes
import mandel_soln
import mandel_gendb

# -------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):

    def setUp(self):
        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": mandel_soln.AnalyticalSoln(),
            "mesh": self.mesh,
        }
        self.checks = [
            Check(
                mesh_entities=["domain"],
                vertex_fields=["displacement"],
                defaults=defaults,
                tolerance=0.5,
            ),
            Check(
                mesh_entities=["domain"],
                vertex_fields=["pressure"],
                defaults=defaults,
                scale=1.0e+6,
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
                vertex_fields = ["displacement"],
                defaults=defaults,
                tolerance=0.5,
            ),
            Check(
                mesh_entities=["poroelastic"],
                vertex_fields = ["pressure"],
                defaults=defaults,
                scale=1.0e+6,
            ),
            Check(
                mesh_entities=["x_neg", "x_pos", "y_neg", "y_pos"],
                filename="output/{name}-{mesh_entity}_info.h5",
                vertex_fields=["initial_amplitude"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["x_neg", "x_pos", "y_neg", "y_pos"],
                vertex_fields=["displacement"],
                defaults=defaults,
                tolerance=0.5,
            ),
            Check(
                mesh_entities=["x_neg", "x_pos", "y_neg", "y_pos"],
                vertex_fields=["pressure"],
                defaults=defaults,
                scale=1.0e+6,
            ),
        ]

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args, mandel_gendb.GenerateDB)


# -------------------------------------------------------------------------------------------------
class TestQuad(TestCase):

    def setUp(self):
        self.name = "mandel_quad"
        self.mesh = meshes.Quad()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["mandel.cfg", "mandel_quad.cfg"])
        return


# -------------------------------------------------------------------------------------------------
class TestTri(TestCase):

    def setUp(self):
        self.name = "mandel_tri"
        self.mesh = meshes.Tri()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["mandel.cfg", "mandel_tri.cfg"])
        return


# ----------------------------------------------------------------------------------------------------------------------
def test_cases():
    return [
        TestQuad,
        TestTri,
    ]


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    FullTestCase.parse_args()

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
