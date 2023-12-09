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
# @file tests/pytests/problems/TestInitialConditionDomain.py
#
# @brief Unit testing of Python InitialConditionDomain object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.problems.InitialConditionDomain import (InitialConditionDomain, initial_conditions)


class TestInitialConditionDomain(TestComponent):
    """Unit testing of InitialConditionDomain object.
    """
    _class = InitialConditionDomain
    _factory = initial_conditions


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestInitialConditionDomain))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
