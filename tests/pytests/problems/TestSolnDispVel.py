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
# @file tests/pytests/problems/TestSolnDispVel.py
#
# @brief Unit testing of Python SolnDispVel object.

import unittest

from pylith.testing.UnitTestApp import (TestAbstractComponent, TestComponent)
from pylith.problems.SolnDispVel import (SolnDispVel, Solution, solution)


class TestSolnDispVel(TestAbstractComponent):
    """Unit testing of SolnDispVel object.
    """
    _class = SolnDispVel


class TestSolutionDispVel(TestAbstractComponent):
    """Unit testing of Solution object.
    """
    _class = Solution
    _factory = solution


if __name__ == "__main__":
    suite = unittest.TestSuite()
    classes = [
        TestSolnDispVel,
        TestSolutionDispVel,
    ]
    for cls in classes:
        suite.addTest(unittest.makeSuite(cls))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
