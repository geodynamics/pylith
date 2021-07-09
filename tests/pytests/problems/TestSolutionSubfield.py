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
# @file tests/pytests/problems/TestSolutionSubfield.py
#
# @brief Unit testing of Python SolutionSubfield object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.problems.SolutionSubfield import (SolutionSubfield, soln_subfield)


class TestSolutionSubfield(TestComponent):
    """Unit testing of SolutionSubfield object.
    """
    _class = SolutionSubfield
    _factory = soln_subfield


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSolutionSubfield))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
