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
