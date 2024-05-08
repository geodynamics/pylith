# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

import unittest

from pylith.testing.TestCases import make_suite

class TestConstants(unittest.TestCase):
    """Unit testing of constants.
    """
  
    def test_maxdouble(self):
        from pylith.utils.utils import maxdouble
        self.assertAlmostEqual(1.0, maxdouble()/1.0e+99, 7)
        return


    def test_maxfloat(self):
        from pylith.utils.utils import maxfloat
        self.assertAlmostEqual(1.0, maxfloat()/1.0e+30, 7)
        return


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestConstants]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file 
