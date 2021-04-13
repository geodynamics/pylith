#!/usr/bin/env nemesis
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#
# @file tests/pytests/meshio/TestDataWriterHDF5.py
#
# @brief Unit testing of Python DataWriterHDF5 object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.meshio.DataWriterHDF5 import (DataWriterHDF5, data_writer)


class TestDataWriterHDF5(TestComponent):
    """Unit testing of DataWriterHDF5 object.
    """
    _class = DataWriterHDF5
    _factory = data_writer


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestDataWriterHDF5))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
