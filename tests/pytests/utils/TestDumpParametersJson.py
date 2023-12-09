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
