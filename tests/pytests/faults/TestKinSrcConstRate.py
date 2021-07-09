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
# @file tests/pytests/faults/TestKinSrcConstRate.py
#
# @brief Unit testing of Python KinSrcConstRate object.

import unittest

from pylith.faults.KinSrcConstRate import KinSrcConstRate
from pylith.tests.UnitTestApp import configureComponent


class TestKinSrcConstRate(unittest.TestCase):
    """Unit testing of KinSrcConstRate object.
    """

    def test_constructor(self):
        src = KinSrcConstRate()

    def test_configure(self):
        src = KinSrcConstRate()
        configureComponent(src)

    def test_factory(self):
        from pylith.faults.KinSrcConstRate import eq_kinematic_src
        src = eq_kinematic_src()
        self.assertTrue(isinstance(src, KinSrcConstRate))


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestKinSrcConstRate))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
