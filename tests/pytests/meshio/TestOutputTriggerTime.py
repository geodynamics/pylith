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
# @file tests/pytests/meshio/TestOutputTriggerTime.py
#
# @brief Unit testing of Python OutputTriggerTime object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.meshio.OutputTriggerTime import (OutputTriggerTime, output_trigger)


class TestOutputTriggerTime(TestComponent):
    """Unit testing of OutputTriggerTime object.
    """
    _class = OutputTriggerTime
    _factory = output_trigger


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestOutputTriggerTime))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
