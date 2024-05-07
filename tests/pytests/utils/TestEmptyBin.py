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

from pylith.testing.TestCases import TestAbstractComponent, make_suite
from pylith.utils.EmptyBin import EmptyBin


class TestEmptyBin(TestAbstractComponent):
    """Unit testing of EmptyBin object.
    """
    _class = EmptyBin


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestEmptyBin]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
