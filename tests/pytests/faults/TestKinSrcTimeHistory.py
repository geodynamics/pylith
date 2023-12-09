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
