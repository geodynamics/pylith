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
# @file tests/pytests/problems/TestInitialConditionPatch.py
#
# @brief Unit testing of Python InitialConditionPatch object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.problems.InitialConditionPatch import (InitialConditionPatch, initial_conditions)


class TestInitialConditionPatch(TestComponent):
    """Unit testing of InitialConditionPatch object.
    """
    _class = InitialConditionPatch
    _factory = initial_conditions


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestInitialConditionPatch))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
