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
# @file tests/pytests/utils/TestNullComponent.py
#
# @brief Unit testing of Python NullComponent object.

import unittest

from pylith.testing.UnitTestApp import TestAbstractComponent
from pylith.utils.NullComponent import NullComponent


class TestNullComponent(TestAbstractComponent):
    """Unit testing of NullComponent object.
    """
    _class = NullComponent


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestNullComponent))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
