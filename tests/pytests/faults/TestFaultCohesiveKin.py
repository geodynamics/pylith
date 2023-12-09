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
# @file tests/pytests/faults/TestFaultCohesiveKin.py
#
# @brief Unit testing of Python FaultCohesiveKin object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.faults.FaultCohesiveKin import (FaultCohesiveKin, fault)


class TestFaultCohesiveKin(TestComponent):
    """Unit testing of FaultCohesiveKin object.
    """
    _class = FaultCohesiveKin
    _factory = fault


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestFaultCohesiveKin))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
