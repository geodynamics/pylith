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

from pylith.testing.TestCases import (TestComponent, configureComponent, make_suite)
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


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestZeroDB]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
