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
