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
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ======================================================================
#
# @file tests/pytests/sources/TestMomentTensorForce.py
#
# @brief Unit testing of Python TestMomentTensorForce object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.sources.Source import (MomentTensorForce, source)


class TestMomentTensorForce(TestComponent):
    """Unit testing of MomentTensorForce object.
    """
    _class = MomentTensorForce
    _factory = source


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestMomentTensorForce))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
