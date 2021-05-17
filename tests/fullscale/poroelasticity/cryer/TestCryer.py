#!/usr/bin/env nemesis
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file tests/fullscale/poroelasticity/cryer/TestCryer.py
#
# @brief Test suite for testing pylith with Cryer's problem.

import unittest

from pylith.tests.FullTestApp import check_data
from pylith.tests.FullTestApp import TestCase as FullTestCase

import meshes
from cryer_soln import AnalyticalSoln
from cryer_gendb import GenerateDB

# We do not include trace_strain in the solution fields, because of the
# poor convergence of the series solution.
SOLUTION_FIELDS = ["displacement", "pressure"]

ratio_tolerance = {'displacement': 1.0, 'pressure': 1.0}
diff_tolerance = {'displacement': 0.5, 'pressure': 0.5}
# ----------------------------------------------------------------------------------------------------------------------


class TestCase(FullTestCase):
    """
    Test suite for testing PyLith with three dimensional poroelasticity
    by means of Cryer's problem.
    """
    DIRICHLET_BOUNDARIES = ["x_neg", "y_neg", "z_neg", "surface_pressure"]
    NEUMANN_BOUNDARIES = ["surface_traction"]
    def setUp(self):
        """
        Setup for test.
        """
        FullTestCase.setUp(self)
        self.exactsoln = AnalyticalSoln()
        return

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args, GenerateDB)
        return

    def test_domain_solution(self):
        filename = "output/{}-domain.h5".format(self.NAME)
        vertexFields = ["displacement", "pressure"]
        check_data(filename, self, self.DOMAIN, vertexFields=vertexFields,
                   ratio_tolerance=ratio_tolerance, diff_tolerance=diff_tolerance)
        return

    def test_material_info(self):
        vertexFields = ['biot_coefficient', 'biot_modulus', 'drained_bulk_modulus',
                        'fluid_density', 'fluid_viscosity', 'isotropic_permeability',
                        'porosity', 'shear_modulus', 'solid_density']
        for material in self.MATERIALS.keys():
            filename = "output/{}-{}_info.h5".format(self.NAME, material)
            check_data(filename, self, self.MATERIALS[material], vertexFields=vertexFields,
                       ratio_tolerance=ratio_tolerance, diff_tolerance=diff_tolerance)
        return

    def test_material_solution(self):
        vertexFields = ["displacement", "pressure"]
        for material in self.MATERIALS.keys():
            filename = "output/{}-{}.h5".format(self.NAME, material)
            check_data(filename, self, self.MATERIALS[material], vertexFields=vertexFields,
                       ratio_tolerance=ratio_tolerance, diff_tolerance=diff_tolerance)
        return

    def test_bcdirichlet_info(self):
        vertexFields = ["initial_amplitude"]
        for bc in self.DIRICHLET_BOUNDARIES:
            self.exactsoln.key = bc
            filename = "output/{}-{}_info.h5".format(self.NAME, bc)
            check_data(filename, self, self.BOUNDARIES[bc], vertexFields=vertexFields,
                       ratio_tolerance=ratio_tolerance, diff_tolerance=diff_tolerance)
        return

    def test_bcdirichlet_solution(self):
        vertexFields = ["displacement", "pressure"]
        for bc in self.DIRICHLET_BOUNDARIES:
            filename = "output/{}-{}.h5".format(self.NAME, bc)
            check_data(filename, self, self.BOUNDARIES[bc], vertexFields=vertexFields,
                       ratio_tolerance=ratio_tolerance, diff_tolerance=diff_tolerance)
        return

# ----------------------------------------------------------------------------------------------------------------------


class TestHex(TestCase, meshes.Hex):
    NAME = "cryer_hex"

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(self, self.NAME, ["cryer.cfg", "cryer_hex.cfg"])
        return


# ----------------------------------------------------------------------------------------------------------------------
class TestTet(TestCase, meshes.Tet):
    NAME = "cryer_tet"

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(self, self.NAME, ["cryer.cfg", "cryer_tet.cfg"])
        return


# ----------------------------------------------------------------------------------------------------------------------
def test_cases():
    return [
        TestHex,
        TestTet,
    ]

# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    FullTestCase.parse_args()

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
