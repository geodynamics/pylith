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
# @file tests/pytests/problems/TestSolnDispPresTracStrain.py
#
# @brief Unit testing of Python SolnDispPresTracStrain object.

import unittest

from pylith.testing.UnitTestApp import (TestAbstractComponent, TestComponent)
from pylith.problems.SolnDispPresTracStrain import (SolnDispPresTracStrain, Solution, solution)


class TestSolnDispPresTracStrain(TestAbstractComponent):
    """Unit testing of SolnDispPresTracStrain object.
    """
    _class = SolnDispPresTracStrain


class TestSolutionDispTracStrain(TestAbstractComponent):
    """Unit testing of Solution object.
    """
    _class = Solution
    _factory = solution


if __name__ == "__main__":
    suite = unittest.TestSuite()
    classes = [
        TestSolnDispPresTracStrain,
        TestSolutionDispTracStrain,
    ]
    for cls in classes:
        suite.addTest(unittest.makeSuite(cls))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
