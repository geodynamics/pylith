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
# @file tests/pytests/meshio/TestOutputTriggerStep.py
#
# @brief Unit testing of Python OutputTriggerStep object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.meshio.OutputTriggerStep import (OutputTriggerStep, output_trigger)


class TestOutputTriggerStep(TestComponent):
    """Unit testing of OutputTriggerStep object.
    """
    _class = OutputTriggerStep
    _factory = output_trigger


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestOutputTriggerStep))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
