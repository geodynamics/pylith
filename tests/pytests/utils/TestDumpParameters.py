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
