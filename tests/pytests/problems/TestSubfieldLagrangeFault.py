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
# @file tests/pytests/problems/TestSubfieldLagrangeFault.py
#
# @brief Unit testing of Python SubfieldLagrangeFault object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.problems.SubfieldLagrangeFault import (SubfieldLagrangeFault, soln_subfield)


class TestSubfieldLagrangeFault(TestComponent):
    """Unit testing of SubfieldLagrangeFault object.
    """
    _class = SubfieldLagrangeFault
    _factory = soln_subfield


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSubfieldLagrangeFault))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
