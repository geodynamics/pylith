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
# @file tests/pytests/problems/TestSubfieldVelocity.py
#
# @brief Unit testing of Python SubfieldVelocity object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.problems.SubfieldVelocity import (SubfieldVelocity, soln_subfield)


class TestSubfieldVelocity(TestComponent):
    """Unit testing of SubfieldVelocity object.
    """
    _class = SubfieldVelocity
    _factory = soln_subfield


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSubfieldVelocity))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
