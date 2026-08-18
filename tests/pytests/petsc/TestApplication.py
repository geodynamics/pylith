#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

import unittest

from pylith.petsc.Application import Application
from pylith.testing.TestCases import make_suite


class TestApplication(unittest.TestCase):
    """Unit testing of Application object.
    """

    def test_constructor(self):
        app = Application()


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestApplication]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
