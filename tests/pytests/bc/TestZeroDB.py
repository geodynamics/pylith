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
# @file tests/pytests/bc/TestZeroDB.py
#
# @brief Unit testing of Python ZeroDB object.

import unittest

from pylith.testing.UnitTestApp import (TestComponent, configureComponent)
from pylith.bc.ZeroDB import (ZeroDB, spatial_database)


class TestZeroDB(TestComponent):
    """Unit testing of ZeroDB object.
    """
    _class = ZeroDB
    _factory = spatial_database

    def test_configure(self):
        zero = ZeroDB()
        configureComponent(zero)
        self.assertEqual(4, len(zero.data))
        for value in zero.data:
            self.assertEqual(0.0, value)
        return


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestZeroDB))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
