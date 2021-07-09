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
# @file tests/pytests/utils/TestDumpParametersJson.py
#
# @brief Unit testing of Python DumpParametersJson object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.utils.DumpParametersJson import (DumpParametersJson, dump_parameters)


class TestDumpParametersJson(TestComponent):
    """Unit testing of DumpParametersJson object.
    """
    _class = DumpParametersJson
    _factory = dump_parameters


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestDumpParametersJson))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
