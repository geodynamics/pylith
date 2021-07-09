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
# @file tests/pytests/meshio/TestOutputSolnPoints.py
#
# @brief Unit testing of Python OutputSolnPoints object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.meshio.OutputSolnPoints import (OutputSolnPoints, observer)


class TestOutputSolnPoints(TestComponent):
    """Unit testing of OutputSolnPoints object.
    """
    _class = OutputSolnPoints
    _factory = observer


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestOutputSolnPoints))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
