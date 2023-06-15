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
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ======================================================================
#
# @file tests/pytests/sources/TestWellboreSource.py
#
# @brief Unit testing of Python TestWellboreSource object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.sources.Source import (WellboreSource, source)


class TestWellboreSource(TestComponent):
    """Unit testing of Elasticity object.
    """
    _class = WellboreSource
    _factory = source


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestWellboreSource))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
