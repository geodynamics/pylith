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
# @file tests/pytests/sources/TestSource.py
#
# @brief Unit testing of Python TestSource object.

import unittest

from pylith.testing.UnitTestApp import TestAbstractComponent
from pylith.sources.Source import Source


class TestSource(TestAbstractComponent):
    """Unit testing of Source object.
    """
    _class = Source


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSource))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
