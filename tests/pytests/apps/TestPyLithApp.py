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
# @file tests/pytests/apps/TestPyLithApp.py
#
# @brief Unit testing of Python PyLithApp object.

import unittest

from pylith.apps.PyLithApp import PyLithApp
from pylith.testing.UnitTestApp import configureComponent


class TestPyLithApp(unittest.TestCase):
    """Unit testing of PyLithApp object.
    """

    def test_constructor(self):
        app = PyLithApp()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestPyLithApp))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
