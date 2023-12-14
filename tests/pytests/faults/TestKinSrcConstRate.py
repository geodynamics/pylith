#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
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
