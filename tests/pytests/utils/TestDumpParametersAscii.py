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
# @file tests/pytests/utils/TestDumpParametersAscii.py
#
# @brief Unit testing of Python DumpParametersAscii object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.utils.DumpParametersAscii import (DumpParametersAscii, dump_parameters)


class TestDumpParametersAscii(TestComponent):
    """Unit testing of DumpParametersAscii object.
    """
    _class = DumpParametersAscii
    _factory = dump_parameters


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestDumpParametersAscii))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
