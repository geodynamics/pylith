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
