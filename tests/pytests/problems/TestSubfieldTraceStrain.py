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
# @file tests/pytests/problems/TestSubfieldTraceStrain.py
#
# @brief Unit testing of Python SubfieldTraceStrain object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.problems.SubfieldTraceStrain import (SubfieldTraceStrain, soln_subfield)


class TestSubfieldTraceStrain(TestComponent):
    """Unit testing of SubfieldTraceStrain object.
    """
    _class = SubfieldTraceStrain
    _factory = soln_subfield


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSubfieldTraceStrain))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
