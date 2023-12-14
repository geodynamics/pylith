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
