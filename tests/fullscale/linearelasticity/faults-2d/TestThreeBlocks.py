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

import unittest

from pylith.testing.FullTestApp import (FullTestCase, Check)

import meshes
import threeblocks_soln


# -------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):

    def setUp(self):
        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": threeblocks_soln.AnalyticalSoln(),
            "mesh": self.mesh,
        }
        self.checks = [
            Check(
                mesh_entities=["domain", "bc_ypos", "points"],
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
                vertex_fields = ["displacement", "cauchy_strain", "cauchy_stress"],
                tolerance = 1.0e-4,
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
        ]

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args)


# -------------------------------------------------------------------------------------------------
class TestQuadGmsh(TestCase):

    def setUp(self):
        self.name = "threeblocks_quad"
        self.mesh = meshes.QuadGmsh()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["threeblocks.cfg", "threeblocks_quad.cfg"])
        return


# -------------------------------------------------------------------------------------------------
class TestTriGmsh(TestCase):

    def setUp(self):
        self.name = "threeblocks_tri"
        self.mesh = meshes.TriGmsh()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["threeblocks.cfg", "threeblocks_tri.cfg"])
        return


# -------------------------------------------------------------------------------------------------
class TestQuadCubit(TestCase):

    def setUp(self):
        self.name = "threeblocks_cubit_quad"
        self.mesh = meshes.QuadCubit()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["threeblocks.cfg", "threeblocks_cubit_quad.cfg"])
        return


# -------------------------------------------------------------------------------------------------
class TestTriCubit(TestCase):

    def setUp(self):
        self.name = "threeblocks_cubit_tri"
        self.mesh = meshes.TriCubit()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["threeblocks.cfg", "threeblocks_cubit_tri.cfg"])
        return


# -------------------------------------------------------------------------------------------------
def test_cases():
    return [
        TestQuadGmsh,
        TestTriGmsh,
        TestQuadCubit,
        TestTriCubit,
    ]


# -------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    FullTestCase.parse_args()

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
