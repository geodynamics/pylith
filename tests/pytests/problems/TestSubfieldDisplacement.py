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
# @file tests/pytests/problems/TestSubfieldDisplacement.py
#
# @brief Unit testing of Python SubfieldDisplacement object.

import unittest

from pylith.testing.TestCases import TestComponent
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
