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

# We cannot test the low-level functionality of the ViscousFriction
# object because it is not exposed to Python. You should really setup
# C++ unit tests using CppUnit as is done for PyLith in addition to
# the simple Python unit tests here.

import unittest


class TestViscousFriction(unittest.TestCase):
    """
    Unit testing of ViscousFriction object.
    """

    def setUp(self):
        """
        Setup test subject.
        """
        from pylith.friction.contrib.ViscousFriction import ViscousFriction
        self.model = ViscousFriction()
        return

    def test_label(self):
        """
        Test constructor.
        """
        label = "viscous friction"
        self.model.label(label)
        self.assertEqual(label, self.model.getLabel())
        return

    def test_timeStep(self):
        """
        Test constructor.
        """
        dt = 2.4
        self.model.timeStep(dt)
        self.assertAlmostEqual(dt, self.model.timeStep(), 5)
        return

    def test_factory(self):
        """
        Test factory method.
        """
        from pylith.friction.contrib.ViscousFriction import friction_model
        f = friction_model()
        return


# End of file
