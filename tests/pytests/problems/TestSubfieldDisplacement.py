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
# @file tests/pytests/problems/TestSubfieldDisplacement.py
#
# @brief Unit testing of Python SubfieldDisplacement object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.problems.SubfieldDisplacement import (SubfieldDisplacement, soln_subfield)


class TestSubfieldDisplacement(TestComponent):
    """Unit testing of SubfieldDisplacement object.
    """
    _class = SubfieldDisplacement
    _factory = soln_subfield


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSubfieldDisplacement))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
