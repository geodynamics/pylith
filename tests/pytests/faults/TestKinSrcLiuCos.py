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
