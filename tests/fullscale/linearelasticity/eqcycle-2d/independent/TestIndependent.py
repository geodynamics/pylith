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
# @file tests/fullscale/linearelasticity/faults-2d/TestTwoBlocks.py
#
# @brief Test suite for testing pylith with 2-D fault shear displacement.

import unittest

from pylith.testing.FullTestApp import check_data
from pylith.testing.FullTestApp import TestCase as FullTestCase

import meshes
from rigid_sliding_soln import AnalyticalSoln as AnalyticalSolnRigidSliding
from simple_shear_soln import AnalyticalSoln as AnalyticalSolnSimpleShear


# ----------------------------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):
    """Test suite for testing PyLith with fault shear displacement.
    """
    FAULTS = ["fault"]
    DIRICHLET_BOUNDARIES = ["bc_xneg", "bc_xpos"]
    NEUMANN_BOUNDARIES = ["bc_yneg", "bc_ypos"]
    OUTPUT_BOUNDARIES = ["bc_ypos"]

    def setUp(self):
        """Setup for test.
        """
        FullTestCase.setUp(self)
        # self.exactsoln = AnalyticalSolnRigidSliding()

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args)

    def test_domain_solution(self):
        filename = "output/{}-domain.h5".format(self.NAME)
        vertexFields = ["displacement"]
        check_data(filename, self, self.DOMAIN, vertexFields=vertexFields)

    def test_material_info(self):
        cellFields = ["density", "bulk_modulus", "shear_modulus"]
        for material in self.MATERIALS:
            filename = "output/{}-{}_info.h5".format(self.NAME, material)
            check_data(filename, self,
                       self.MATERIALS[material], cellFields=cellFields)

    def test_material_solution(self):
        vertexFields = ["displacement"]
        for material in self.MATERIALS:
            filename = "output/{}-{}.h5".format(self.NAME, material)
            check_data(filename, self,
                       self.MATERIALS[material], vertexFields=vertexFields)

    def test_bcdirichlet_info(self):
        vertexFields = ["initial_amplitude"]
        for bc in self.DIRICHLET_BOUNDARIES:
            filename = "output/{}-{}_info.h5".format(self.NAME, bc)
            check_data(filename, self,
                       self.BOUNDARIES[bc], vertexFields=vertexFields)

    def test_bcdirichlet_solution(self):
        vertexFields = ["displacement"]
        for bc in self.DIRICHLET_BOUNDARIES:
            filename = "output/{}-{}.h5".format(self.NAME, bc)
            check_data(filename, self,
                       self.BOUNDARIES[bc], vertexFields=vertexFields)

    def test_boundary_solution(self):
        vertexFields = ["displacement"]
        if "simple_shear" in self.NAME:
            return
        else:
            for bc in self.OUTPUT_BOUNDARIES:
                filename = "output/{}-{}.h5".format(self.NAME, bc)
                check_data(filename, self,
                           self.BOUNDARIES[bc], vertexFields=vertexFields)

    def test_bcneumann_info(self):
        vertexFields = ["initial_amplitude"]
        if "rigid_sliding" in self.NAME:
            return 
        else:
            for bc in self.NEUMANN_BOUNDARIES:
                filename = "output/{}-{}_info.h5".format(self.NAME, bc)
                check_data(filename, self,
                           self.BOUNDARIES[bc], vertexFields=vertexFields)

    def test_bcneumann_solution(self):
        vertexFields = ["displacement"]
        if "rigid_sliding" in self.NAME:
            return
        else:
            for bc in self.NEUMANN_BOUNDARIES:
                filename = "output/{}-{}.h5".format(self.NAME, bc)
                check_data(filename, self,
                           self.BOUNDARIES[bc], vertexFields=vertexFields)

# ----------------------------------------------------------------------------------------------------------------------
class RigidSliding_Quad(TestCase, meshes.Quad):
    NAME = "rigid_sliding_quad"

    def setUp(self):
        TestCase.setUp(self)
        self.exactsoln = AnalyticalSolnRigidSliding()
        TestCase.run_pylith(
            self, self.NAME, ["independent.cfg", "independent_quad.cfg"])
        return


# ----------------------------------------------------------------------------------------------------------------------
class RigidSliding_Tri(TestCase, meshes.Tri):
    NAME = "rigid_sliding_tri"

    def setUp(self):
        TestCase.setUp(self)
        self.exactsoln = AnalyticalSolnRigidSliding()
        TestCase.run_pylith(
            self, self.NAME, ["independent.cfg", "independent_tri.cfg"])
        return

# ----------------------------------------------------------------------------------------------------------------------
class SimpleShear_Quad(TestCase, meshes.Quad):
    NAME = "simple_shear_quad"

    def setUp(self):
        TestCase.setUp(self)
        self.exactsoln = AnalyticalSolnSimpleShear()
        return


# ----------------------------------------------------------------------------------------------------------------------
class SimpleShear_Tri(TestCase, meshes.Tri):
    NAME = "simple_shear_tri"

    def setUp(self):
        TestCase.setUp(self)
        self.exactsoln = AnalyticalSolnSimpleShear()
        return

# ----------------------------------------------------------------------------------------------------------------------
def test_cases():
    return [
        RigidSliding_Quad,
        SimpleShear_Quad,
        RigidSliding_Tri,
        SimpleShear_Tri
    ]


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    FullTestCase.parse_args()

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
