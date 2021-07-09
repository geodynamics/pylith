#!/usr/bin/env nemesis
#
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
# @file tests/pytests/faults/TestKinSrcRamp.py
#
# @brief Unit testing of Python KinSrcRamp object.

import unittest

from pylith.faults.KinSrcRamp import KinSrcRamp
from pylith.tests.UnitTestApp import configureComponent


class TestKinSrcRamp(unittest.TestCase):
    """Unit testing of KinSrcRamp object.
    """

    def test_constructor(self):
        src = KinSrcRamp()

    def test_configure(self):
        src = KinSrcRamp()
        configureComponent(src)

    def test_factory(self):
        from pylith.faults.KinSrcRamp import eq_kinematic_src
        src = eq_kinematic_src()
        self.assertTrue(isinstance(src, KinSrcRamp))


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestKinSrcRamp))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
