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
