#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
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
