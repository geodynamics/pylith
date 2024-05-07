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
# @file tests/pytests/problems/TestSubfieldTraceStrain.py
#
# @brief Unit testing of Python SubfieldTraceStrain object.

import unittest

from pylith.testing.TestCases import TestComponent
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
