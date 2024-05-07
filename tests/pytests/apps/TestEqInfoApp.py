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

import unittest

from pylith.apps.EqInfoApp import EqInfoApp
from pylith.testing.TestCases import make_suite


class TestEqInfoApp(unittest.TestCase):
    """Unit testing of EqInfoApp object.
    """

    def test_constructor(self):
        app = EqInfoApp()


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestEqInfoApp]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
