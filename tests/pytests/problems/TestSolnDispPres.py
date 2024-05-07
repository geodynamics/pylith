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
# @file tests/pytests/problems/TestSolnDispPres.py
#
# @brief Unit testing of Python SolnDispPres object.

import unittest

from pylith.testing.TestCases import (TestAbstractComponent, TestComponent)
from pylith.problems.SolnDispPres import (SolnDispPres, Solution, solution)


class TestSolnDispPres(TestAbstractComponent):
    """Unit testing of SolnDispPres object.
    """
    _class = SolnDispPres


class TestSolutionDispPres(TestAbstractComponent):
    """Unit testing of Solution object.
    """
    _class = Solution
    _factory = solution


if __name__ == "__main__":
    suite = unittest.TestSuite()
    classes = [
        TestSolnDispPres,
        TestSolutionDispPres,
    ]
    for cls in classes:
        suite.addTest(unittest.makeSuite(cls))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
