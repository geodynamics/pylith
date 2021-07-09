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
