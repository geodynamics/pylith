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
# @file tests/pytests/problems/TestSolnDispLagrange.py
#
# @brief Unit testing of Python SolnDispLagrange object.

import unittest

from pylith.testing.UnitTestApp import (TestAbstractComponent, TestComponent)
from pylith.problems.SolnDispLagrange import (SolnDispLagrange, Solution, solution)


class TestSolnDispLagrange(TestAbstractComponent):
    """Unit testing of SolnDispLagrange object.
    """
    _class = SolnDispLagrange


class TestSolutionDispLagrange(TestAbstractComponent):
    """Unit testing of Solution object.
    """
    _class = Solution
    _factory = solution


if __name__ == "__main__":
    suite = unittest.TestSuite()
    classes = [
        TestSolnDispLagrange,
        TestSolutionDispLagrange,
    ]
    for cls in classes:
        suite.addTest(unittest.makeSuite(cls))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
