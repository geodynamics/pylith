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
# @file tests/pytests/faults/TestFaultCohesiveImpulses.py
#
# @brief Unit testing of Python FaultCohesiveImpulses object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.faults.FaultCohesiveImpulses import (FaultCohesiveImpulses, fault)


class TestFaultCohesiveImpulses(TestComponent):
    """Unit testing of FaultCohesiveImpulses object.
    """
    _class = FaultCohesiveImpulses
    _factory = fault


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestFaultCohesiveImpulses))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
