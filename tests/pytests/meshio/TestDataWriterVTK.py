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
