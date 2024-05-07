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
# @file tests/pytests/problems/TestSolnDisp.py
#
# @brief Unit testing of Python SolnDisp object.

import unittest

from pylith.testing.TestCases import TestAbstractComponent
from pylith.problems.SolnDisp import SolnDisp


class TestSolnDisp(TestAbstractComponent):
    """Unit testing of SolnDisp object.
    """
    _class = SolnDisp


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSolnDisp))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
