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
# @file tests/pytests/utils/TestDumpParameters.py
#
# @brief Unit testing of Python DumpParameters object.

import unittest

from pylith.testing.UnitTestApp import TestAbstractComponent
from pylith.utils.DumpParameters import DumpParameters


class TestDumpParameters(TestAbstractComponent):
    """Unit testing of DumpParameters object.
    """
    _class = DumpParameters


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestDumpParameters))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
