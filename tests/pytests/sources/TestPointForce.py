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
# @file tests/pytests/sources/TestPointForce.py
#
# @brief Unit testing of Python TestPointForce object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.sources.Source import (PointForce, source)


class TestPointForce(TestComponent):
    """Unit testing of Elasticity object.
    """
    _class = PointForce
    _factory = source


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestPointForce))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
