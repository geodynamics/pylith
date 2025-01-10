# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

import unittest

from pylith.testing.TestCases import TestComponent, make_suite
from pylith.meshio.DataWriterHDF5Ext import (DataWriterHDF5Ext, data_writer)


class TestDataWriterHDF5Ext(TestComponent):
    """Unit testing of DataWriterHDF5Ext object.
    """
    _class = DataWriterHDF5Ext
    _factory = data_writer


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestDataWriterHDF5Ext]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
