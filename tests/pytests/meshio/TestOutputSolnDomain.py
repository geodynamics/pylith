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
# @file tests/pytests/meshio/TestOutputSolnDomain.py
#
# @brief Unit testing of Python OutputSolnDomain object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.meshio.OutputSolnDomain import (OutputSolnDomain, observer)


class TestOutputSolnDomain(TestComponent):
    """Unit testing of OutputSolnDomain object.
    """
    _class = OutputSolnDomain
    _factory = observer


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestOutputSolnDomain))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
