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
# @file tests/pytests/problems/TestSolution.py
#
# @brief Unit testing of Python Solution object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.problems.Solution import (Solution, solution)


class TestSolution(TestComponent):
    """Unit testing of Solution object.
    """
    _class = Solution
    _factory = solution


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSolution))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
