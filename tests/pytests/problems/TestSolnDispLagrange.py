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
# @file tests/pytests/problems/TestSolnDispLagrange.py
#
# @brief Unit testing of Python SolnDispLagrange object.

import unittest

from pylith.testing.TestCases import (TestAbstractComponent, TestComponent)
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
