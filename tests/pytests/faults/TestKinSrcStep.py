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
# @file tests/pytests/faults/TestKinSrcStep.py
#
# @brief Unit testing of Python KinSrcStep object.

import unittest

from pylith.faults.KinSrcStep import KinSrcStep
from pylith.tests.UnitTestApp import configureComponent


class TestKinSrcStep(unittest.TestCase):
    """Unit testing of KinSrcStep object.
    """

    def test_constructor(self):
        src = KinSrcStep()

    def test_configure(self):
        src = KinSrcStep()
        configureComponent(src)

    def test_factory(self):
        from pylith.faults.KinSrcStep import eq_kinematic_src
        src = eq_kinematic_src()
        self.assertTrue(isinstance(src, KinSrcStep))


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestKinSrcStep))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
