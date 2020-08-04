#!/usr/bin/env python
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
# @file tests/fullscale/linearporoelasticity/terzaghi/TestTerzaghi.py
#
# @brief Test suite for testing pylith with Terzaghi's problem.

import unittest

from pylith.tests.FullTestApp import check_data
from pylith.tests.FullTestApp import TestCase as FullTestCase

import meshes
from terzaghi_soln import AnalyticalSoln
from terzaghi_gendb import GenerateDB

tolerance = 1e-2
# ----------------------------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):
    """
    Test suite for testing PyLith with one dimensional poroelasticity
    by means of Terzaghi's problem.
    """
    DIRICHLET_BOUNDARIES = ["x_neg", "x_pos", "y_neg_dir", "y_pos"]
    NEUMANN_BOUNDARIES = ["y_neg_neu"]

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
        vertexFields = ["displacement", "pressure", "trace_strain"]
        check_data(filename, self, self.DOMAIN, vertexFields=vertexFields, tolerance=tolerance)
        return

    def test_material_info(self):
        cellFields = ["solid_density", "fluid_density", "fluid_viscosity", "shear_modulus", "undrained_bulk_modulus", "biot_coefficient", "biot_modulus", "isotropic_permeability"]
        for material in self.MATERIALS.keys():
            filename = "output/{}-{}_info.h5".format(self.NAME, material)
            check_data(filename, self, self.MATERIALS[material], cellFields=cellFields, tolerance=tolerance)
        return

    def test_material_solution(self):
        vertexFields = ["displacement", "pressure", "trace_strain"]
        for material in self.MATERIALS.keys():
            filename = "output/{}-{}.h5".format(self.NAME, material)
            check_data(filename, self, self.MATERIALS[material], vertexFields=vertexFields, tolerance=tolerance)
        return

#    def test_bcdirichlet_info(self):
#        vertexFields = ["initial_amplitude"]
#        for bc in self.DIRICHLET_BOUNDARIES:
#            self.exactsoln.key = bc
#            filename = "output/{}-{}_info.h5".format(self.NAME, bc)
#            check_data(filename, self, self.BOUNDARIES[bc], vertexFields=vertexFields)
#        return

#    def test_bcdirichlet_solution(self):
#        vertexFields = ["displacement", "pressure", "trace_strain"]
#        for bc in self.DIRICHLET_BOUNDARIES:
#            filename = "output/{}-{}.h5".format(self.NAME, bc)
#            check_data(filename, self, self.BOUNDARIES[bc], vertexFields=vertexFields)
#        return

#    def test_bcneumann_info(self):
#        vertexFields = ["initial_amplitude"]
#        for bc in self.NEUMANN_BOUNDARIES:
#            self.exactsoln.key = bc
#            filename = "output/{}-{}_info.h5".format(self.NAME, bc)
#            check_data(filename, self, self.BOUNDARIES[bc], vertexFields=vertexFields)
#        return

#    def test_bcneumann_solution(self):
#        vertexFields = ["displacement", "pressure", "trace_strain"]
#        for bc in self.NEUMANN_BOUNDARIES:
#            filename = "output/{}-{}.h5".format(self.NAME, bc)
#            check_data(filename, self, self.BOUNDARIES[bc], vertexFields=vertexFields)
#        return

# ----------------------------------------------------------------------------------------------------------------------
class TestQuad(TestCase, meshes.Quad):
    NAME = "terzaghi_quad"

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(self, self.NAME, ["terzaghi.cfg", "terzaghi_quad.cfg"])
        return


# ----------------------------------------------------------------------------------------------------------------------
class TestTri(TestCase, meshes.Tri):
    NAME = "terzaghi_tri"

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(self, self.NAME, ["terzaghi.cfg", "terzaghi_tri.cfg"])
        return


# ----------------------------------------------------------------------------------------------------------------------
def test_cases():
    return [
        TestQuad,
#        TestTri,
    ]


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    FullTestCase.parse_args()

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
