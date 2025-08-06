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

from pylith.testing.FullTestApp import FullTestCase, Check
from pylith.testing import TestCases

import meshes
import gravity_soln


# -------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):

    def setUp(self):
        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": gravity_soln.AnalyticalSoln(),
            "mesh": self.mesh,
        }
        self.checks = [
            Check(
                mesh_entities=["domain"],
                vertex_fields=["displacement", "pressure", "trace_strain"],
                final_time_only=True,
                defaults=defaults,
            ),
            Check(
                mesh_entities=["poroelastic"],
                filename="output/{name}-{mesh_entity}_info.h5",
                cell_fields=[
                    "solid_density",
                    "fluid_density",
                    "fluid_viscosity",
                    "porosity",
                    "shear_modulus",
                    "drained_bulk_modulus",
                    "biot_coefficient",
                    "biot_modulus",
                    "isotropic_permeability",
                ],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["poroelastic"],
                vertex_fields=[
                    "displacement",
                    "pressure",
                    "trace_strain",
                    "cauchy_strain",
                    "cauchy_stress",
                ],
                final_time_only=True,
                defaults=defaults,
            ),
            Check(
                mesh_entities=[
                    "bc_disp_xneg",
                    "bc_disp_xpos",
                    "bc_disp_yneg",
                ],
                filename="output/{name}-{mesh_entity}_info.h5",
                vertex_fields=["initial_amplitude", "normal_dir", "tangential_dir"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=[
                    "bc_disp_xneg",
                    "bc_disp_xpos",
                    "bc_disp_yneg",
                    "bc_press_ypos",
                ],
                vertex_fields=["displacement", "pressure", "trace_strain"],
                final_time_only=True,
                defaults=defaults,
            ),
        ]

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args)


# -------------------------------------------------------------------------------------------------
class TestQuadGmsh(TestCase):

    def setUp(self):
        self.name = "gravity_quad"
        self.mesh = meshes.QuadGmsh()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["gravity.cfg", "gravity_quad.cfg"])
        return


# -------------------------------------------------------------------------------------------------
class TestTriGmsh(TestCase):

    def setUp(self):
        self.name = "gravity_tri"
        self.mesh = meshes.TriGmsh()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["gravity.cfg", "gravity_tri.cfg"])
        return


# -------------------------------------------------------------------------------------------------
def load_tests(loader, tests, pattern):
    TEST_CLASSES = (
        TestQuadGmsh,
        TestTriGmsh,
    )
    return TestCases.make_suite(test_classes=TEST_CLASSES, loader=loader)


# -------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    FullTestCase.parse_args()
    unittest.main(verbosity=2, argv=["test"])


# End of file
