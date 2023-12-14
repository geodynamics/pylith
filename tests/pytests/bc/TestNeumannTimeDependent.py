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
# @file tests/pytests/bc/TestNeumannTimeDependent.py
#
# @brief Unit testing of Python NeumannTimeDependent object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.bc.NeumannTimeDependent import (NeumannTimeDependent, boundary_condition)


class TestNeumannTimeDependent(TestComponent):
    """Unit testing of NeumannTimeDependent object.
    """
    _class = NeumannTimeDependent
    _factory = boundary_condition


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestNeumannTimeDependent))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
