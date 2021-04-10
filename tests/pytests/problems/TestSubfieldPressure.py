#!/usr/bin/env nemesis
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#
# @file tests/pytests/problems/TestSubfieldPressure.py
#
# @brief Unit testing of Python SubfieldPressure object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.problems.SubfieldPressure import (SubfieldPressure, soln_subfield)


class TestSubfieldPressure(TestComponent):
    """Unit testing of SubfieldPressure object.
    """
    _class = SubfieldPressure
    _factory = soln_subfield


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSubfieldPressure))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
