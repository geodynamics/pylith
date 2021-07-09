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
