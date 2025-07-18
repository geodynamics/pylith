#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

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
                mesh_entities=["bc_xneg", "bc_xpos", "bc_yneg", "bc_ypos"],
                filename="output/{name}-{mesh_entity}_info.h5",
                vertex_fields=["initial_amplitude", "normal_dir", "tangential_dir"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_xpos", "bc_yneg", "bc_ypos"],
                vertex_fields=["displacement"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["fault_xneg", "fault_xmid"],
                filename="output/{name}-{mesh_entity}_info.h5",
                vertex_fields = ["normal_dir", "strike_dir"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["fault_xneg", "fault_xmid"],
                vertex_fields=["slip", "traction_change"],
                defaults=defaults,
            ),
        ]

    def run_pylith(self, testName, args, nprocs=1):
        FullTestCase.run_pylith(self, testName, args, shearnoslip_gendb.GenerateDB, nprocs=nprocs)


# -------------------------------------------------------------------------------------------------
class TestQuadGmsh(TestCase):

    def setUp(self):
        self.name = "shearnoslip_quad"
        self.mesh = meshes.QuadGmsh()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["shearnoslip.cfg", "shearnoslip_quad.cfg"], nprocs=1)


# -------------------------------------------------------------------------------------------------
class TestTriGmsh(TestCase):

    def setUp(self):
        self.name = "shearnoslip_tri"
        self.mesh = meshes.TriGmsh()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["shearnoslip.cfg", "shearnoslip_tri.cfg"], nprocs=1)


# -------------------------------------------------------------------------------------------------
class TestQuadGmshRefineOutput(TestCase):

    def setUp(self):
        self.name = "shearnoslip_refineoutput_quad"
        self.mesh = meshes.QuadGmshRefineOutput()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["shearnoslip.cfg", "shearnoslip_refineoutput.cfg", "shearnoslip_refineoutput_quad.cfg"], nprocs=2)


# -------------------------------------------------------------------------------------------------
class TestTriGmshRefineOutput(TestCase):

    def setUp(self):
        self.name = "shearnoslip_refineoutput_tri"
        self.mesh = meshes.TriGmshRefineOutput()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["shearnoslip.cfg", "shearnoslip_refineoutput.cfg", "shearnoslip_refineoutput_tri.cfg"], nprocs=1)


# -------------------------------------------------------------------------------------------------
def test_cases():
    return [
        TestQuadGmsh,
        TestTriGmsh,
        TestQuadGmshRefineOutput,
        TestTriGmshRefineOutput,
    ]


# -------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    FullTestCase.parse_args()

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
