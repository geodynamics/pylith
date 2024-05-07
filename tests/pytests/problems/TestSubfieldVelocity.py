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
# @file tests/pytests/problems/TestSubfieldVelocity.py
#
# @brief Unit testing of Python SubfieldVelocity object.

import unittest

from pylith.testing.TestCases import TestComponent
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
