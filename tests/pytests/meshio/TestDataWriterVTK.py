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
# @file tests/pytests/meshio/TestDataWriterVTK.py
#
# @brief Unit testing of Python DataWriterVTK object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.meshio.DataWriterVTK import (DataWriterVTK, data_writer)


class TestDataWriterVTK(TestComponent):
    """Unit testing of DataWriterVTK object.
    """
    _class = DataWriterVTK
    _factory = data_writer


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestDataWriterVTK))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
