# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

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
