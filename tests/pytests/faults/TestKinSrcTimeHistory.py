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
# @file tests/pytests/faults/TestKinSrcTimeHistory.py
#
# @brief Unit testing of Python KinSrcTimeHistory object.

import unittest

from pylith.faults.KinSrcTimeHistory import KinSrcTimeHistory
from pylith.tests.UnitTestApp import configureComponent


class TestKinSrcTimeHistory(unittest.TestCase):
    """Unit testing of KinSrcTimeHistory object.
    """

    def test_constructor(self):
        src = KinSrcTimeHistory()

    def test_configure(self):
        src = KinSrcTimeHistory()
        configureComponent(src)

    def test_factory(self):
        from pylith.faults.KinSrcTimeHistory import eq_kinematic_src
        src = eq_kinematic_src()
        self.assertTrue(isinstance(src, KinSrcTimeHistory))


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestKinSrcTimeHistory))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
