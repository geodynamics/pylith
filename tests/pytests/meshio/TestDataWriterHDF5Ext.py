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
# @file tests/pytests/meshio/TestDataWriterHDF5Ext.py
#
# @brief Unit testing of Python DataWriterHDF5Ext object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.meshio.DataWriterHDF5Ext import (DataWriterHDF5Ext, data_writer)


class TestDataWriterHDF5Ext(TestComponent):
    """Unit testing of DataWriterHDF5Ext object.
    """
    _class = DataWriterHDF5Ext
    _factory = data_writer


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestDataWriterHDF5Ext))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
