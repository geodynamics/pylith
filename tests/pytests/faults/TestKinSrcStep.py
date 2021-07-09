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
