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
# @file tests/pytests/faults/TestKinSrcBrune.py
#
# @brief Unit testing of Python KinSrcBrune object.

import unittest

from pylith.faults.KinSrcBrune import KinSrcBrune
from pylith.tests.UnitTestApp import configureComponent


class TestKinSrcBrune(unittest.TestCase):
    """Unit testing of KinSrcBrune object.
    """

    def test_constructor(self):
        src = KinSrcBrune()

    def test_configure(self):
        src = KinSrcBrune()
        configureComponent(src)

    def test_factory(self):
        from pylith.faults.KinSrcBrune import eq_kinematic_src
        src = eq_kinematic_src()
        self.assertTrue(isinstance(src, KinSrcBrune))


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestKinSrcBrune))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
