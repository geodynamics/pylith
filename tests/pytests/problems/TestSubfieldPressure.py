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
