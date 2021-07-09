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
# @file tests/pytests/problems/TestTimeDependent.py
#
# @brief Unit testing of Python TimeDependent object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.problems.TimeDependent import (TimeDependent, problem)


class TestTimeDependent(TestComponent):
    """Unit testing of TimeDependent object.
    """
    _class = TimeDependent
    _factory = problem


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTimeDependent))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
