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
# @file tests/pytests/apps/TestPyLithApp.py
#
# @brief Unit testing of Python PyLithApp object.

import unittest

from pylith.apps.PyLithApp import PyLithApp
from pylith.testing.TestCases import make_suite


class TestPyLithApp(unittest.TestCase):
    """Unit testing of PyLithApp object.
    """

    def test_constructor(self):
        app = PyLithApp()


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestPyLithApp]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
