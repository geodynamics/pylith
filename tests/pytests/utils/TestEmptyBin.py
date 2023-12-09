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
# @file tests/pytests/utils/TestEmptyBin.py
#
# @brief Unit testing of Python EmptyBin object.

import unittest

from pylith.testing.UnitTestApp import TestAbstractComponent
from pylith.utils.EmptyBin import EmptyBin


class TestEmptyBin(TestAbstractComponent):
    """Unit testing of EmptyBin object.
    """
    _class = EmptyBin


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestEmptyBin))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
