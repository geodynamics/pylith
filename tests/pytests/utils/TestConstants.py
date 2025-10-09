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

from pylith.testing.TestCases import make_suite


class TestConstants(unittest.TestCase):
    """Unit testing of constants."""

    def test_g_acc(self):
        from pylith.utils.utils import g_acc

        self.assertAlmostEqual(1.0, g_acc() / 9.80665, 7)

    def test_max_double(self):
        from pylith.utils.utils import max_double

        assert max_double() > 1.0e99

    def test_max_float(self):
        from pylith.utils.utils import max_float

        assert max_float() > 1.0e30


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestConstants]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
