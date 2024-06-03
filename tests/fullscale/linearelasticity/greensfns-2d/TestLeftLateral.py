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

import unittest

from pylith.testing.FullTestApp import (FullTestCase, Check, check_sizes)

import meshes
import faultimpulses_soln

# -------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):

    def setUp(self):
        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": faultimpulses_soln.AnalyticalSoln(),
            "mesh": self.mesh,
        }
        self.checks = [
            Check(
                mesh_entities=["mat_xneg", "mat_xpos"],
                filename="output/{name}-{mesh_entity}_info.h5",
                cell_fields = ["density", "bulk_modulus", "shear_modulus"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_xpos"],
                filename="output/{name}-{mesh_entity}_info.h5",
                vertex_fields=["initial_amplitude"],
                defaults=defaults,
            ),
        ]

        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": faultimpulses_soln.SolnDims(self.nimpulses, self.basis_order),
            "mesh": self.mesh,
        }
        self.dim_checks = [
            Check(
                mesh_entities=["domain", "bc_ypos"],
                vertex_fields=["displacement"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["mat_xneg", "mat_xpos"],
                vertex_fields = ["displacement"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_xpos"],
                vertex_fields=["displacement"],
                defaults=defaults,
            ),
        ]
        if self.basis_order["slip"] == 0:
            self.dim_checks.append(
                Check(
                    mesh_entities=["fault"],
                    cell_fields=["slip"],
                    defaults=defaults,
                )
            )
        else:
            self.dim_checks.append(
                Check(
                    mesh_entities=["fault"],
                    vertex_fields=["slip"],
                    defaults=defaults,
                )
            )

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args)

    def test_output_size(self):
        for check in self.dim_checks:
            for mesh_entity in check.mesh_entities:
                filename = check.filename.format(name=self.name, mesh_entity=mesh_entity)
                with self.subTest(filename=filename):
                    check_sizes(self, filename, check, mesh_entity, check.mesh.ENTITIES[mesh_entity])



# -------------------------------------------------------------------------------------------------
class TestQuadBasis0(TestCase):

    def setUp(self):
        self.name = "leftlateral_b0_quad"
        self.mesh = meshes.QuadGmsh()
        self.nimpulses = 9
        self.basis_order = {
            "displacement": 1,
            "slip": 0
        }
        super().setUp()

        TestCase.run_pylith(self, self.name, ["leftlateral_b0.cfg", "leftlateral_b0_quad.cfg"])


# -------------------------------------------------------------------------------------------------
class TestQuadBasis1(TestCase):

    def setUp(self):
        self.name = "leftlateral_b1_quad"
        self.mesh = meshes.QuadGmsh()
        self.nimpulses = 10
        self.basis_order = {
            "displacement": 1,
            "slip": 1,
        }
        super().setUp()

        TestCase.run_pylith(self, self.name, ["leftlateral_b1.cfg", "leftlateral_b1_quad.cfg"])


# -------------------------------------------------------------------------------------------------
class TestQuadBasis2(TestCase):

    def setUp(self):
        self.name = "leftlateral_b2_quad"
        self.mesh = meshes.QuadGmsh()
        self.nimpulses = 19
        self.basis_order = {
            "displacement": 1,
            "slip": 1,
        }
        super().setUp()

        TestCase.run_pylith(self, self.name, ["leftlateral_b2.cfg", "leftlateral_b2_quad.cfg"])


# -------------------------------------------------------------------------------------------------
class TestTriBasis0(TestCase):

    def setUp(self):
        self.name = "leftlateral_b0_tri"
        self.mesh = meshes.TriGmsh()
        self.nimpulses = 8
        self.basis_order = {
            "displacement": 1,
            "slip": 0
        }
        super().setUp()

        TestCase.run_pylith(self, self.name, ["leftlateral_b0.cfg", "leftlateral_b0_tri.cfg"])


# -------------------------------------------------------------------------------------------------
class TestTriBasis1(TestCase):

    def setUp(self):
        self.name = "leftlateral_b1_tri"
        self.mesh = meshes.TriGmsh()
        self.nimpulses = 9
        self.basis_order = {
            "displacement": 1,
            "slip": 1,
        }
        super().setUp()

        TestCase.run_pylith(self, self.name, ["leftlateral_b1.cfg", "leftlateral_b1_tri.cfg"])


# -------------------------------------------------------------------------------------------------
class TestTriBasis2(TestCase):

    def setUp(self):
        self.name = "leftlateral_b2_tri"
        self.mesh = meshes.TriGmsh()
        self.nimpulses = 17
        self.basis_order = {
            "displacement": 1,
            "slip": 1,
        }
        super().setUp()

        TestCase.run_pylith(self, self.name, ["leftlateral_b2.cfg", "leftlateral_b2_tri.cfg"])


# -------------------------------------------------------------------------------------------------
def test_cases():
    return [
        TestQuadBasis0,
        TestQuadBasis1,
        TestQuadBasis2,
        TestTriBasis0,
        TestTriBasis1,
        TestTriBasis2,
    ]


# -------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    FullTestCase.parse_args()

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
