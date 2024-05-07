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
# @file tests/pytests/problems/TestSolnDispVelLagrange.py
#
# @brief Unit testing of Python SolnDispVelLagrange object.

import unittest

from pylith.testing.TestCases import (TestAbstractComponent, TestComponent)
from pylith.problems.SolnDispVelLagrange import (SolnDispVelLagrange, Solution, solution)


class TestSolnDispVelLagrange(TestAbstractComponent):
    """Unit testing of SolnDispVelLagrange object.
    """
    _class = SolnDispVelLagrange


class TestSolutionDispVelLagrange(TestAbstractComponent):
    """Unit testing of Solution object.
    """
    _class = Solution
    _factory = solution


if __name__ == "__main__":
    suite = unittest.TestSuite()
    classes = [
        TestSolnDispVelLagrange,
        TestSolutionDispVelLagrange,
    ]
    for cls in classes:
        suite.addTest(unittest.makeSuite(cls))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
