#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

import unittest

from pylith.testing.FullTestApp import TestDriver, FullTestCase

import TestBodyForce
import TestGravity
import TestGravityBodyForce

TEST_MODULES = (
    TestBodyForce,
    TestGravity,
    TestGravityBodyForce,
)


class TestApp(TestDriver):
    """Driver application for full-scale tests."""

    def __init__(self):
        """Constructor."""
        TestDriver.__init__(self)

    def _suite(self):
        """Create test suite."""

        loader = unittest.defaultTestLoader
        suite = unittest.TestSuite()
        for mod in TEST_MODULES:
            suite.addTests(loader.loadTestsFromModule(mod))
        return suite


# ----------------------------------------------------------------------
if __name__ == "__main__":
    FullTestCase.parse_args()
    TestApp().main()


# End of file
