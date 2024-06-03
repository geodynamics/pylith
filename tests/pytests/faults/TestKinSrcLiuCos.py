#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file tests/pytests/faults/TestKinSrcLiuCos.py
#
# @brief Unit testing of Python KinSrcLiuCos object.

import unittest

from pylith.faults.KinSrcLiuCos import KinSrcLiuCos
from pylith.tests.UnitTestApp import configureComponent


class TestKinSrcLiuCos(unittest.TestCase):
    """Unit testing of KinSrcLiuCos object.
    """

    def test_constructor(self):
        src = KinSrcLiuCos()

    def test_configure(self):
        src = KinSrcLiuCos()
        configureComponent(src)

    def test_factory(self):
        from pylith.faults.KinSrcLiuCos import eq_kinematic_src
        src = eq_kinematic_src()
        self.assertTrue(isinstance(src, KinSrcLiuCos))


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestKinSrcLiuCos))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
