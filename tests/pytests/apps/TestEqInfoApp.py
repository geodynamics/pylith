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
# @file tests/pytests/apps/TestEqInfoApp.py
#
# @brief Unit testing of Python EqInfoApp object.

import unittest

from pylith.apps.EqInfoApp import EqInfoApp
from pylith.testing.UnitTestApp import configureComponent


class TestEqInfoApp(unittest.TestCase):
    """Unit testing of EqInfoApp object.
    """

    def test_constructor(self):
        app = EqInfoApp()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestEqInfoApp))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
