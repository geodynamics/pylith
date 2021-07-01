# ======================================================================
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
# ======================================================================
#

# We cannot test the low-level functionality of the PlaneStrainState
# object because it is not exposed to Python. You should really setup
# C++ unit tests using CppUnit as is done for PyLith in addition to
# the simple Python unit tests here.

import unittest


class TestPlaneStrainState(unittest.TestCase):
    """
    Unit testing of PlaneStrainState object.
    """

    def setUp(self):
        """
        Setup test subject.
        """
        from pylith.materials.contrib.PlaneStrainState import PlaneStrainState
        self.material = PlaneStrainState()
        return

    def test_constructor(self):
        """
        Test constructor.
        """
        self.assertEqual(2, self.material.dimension())
        return

    def test_useElasticBehavior(self):
        """
        Test useElasticBehavior().
        """
        self.material.useElasticBehavior(False)
        return

    def testHasStateVars(self):
        self.failUnless(self.material.hasStateVars())
        return

    def testTensorSize(self):
        self.assertEqual(3, self.material.tensorSize())
        return

    def testNeedNewJacobian(self):
        """
        Test needNewJacobian().
        """
        # Default should be False.
        self.failIf(self.material.needNewJacobian())

        # Changing time step should not require new Jacobian.
        self.material.timeStep(1.0)
        self.material.timeStep(2.0)
        self.failIf(self.material.needNewJacobian())
        return

    def test_factory(self):
        """
        Test factory method.
        """
        from pylith.materials.contrib.PlaneStrainState import material
        m = material()
        return


# End of file
