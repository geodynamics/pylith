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

from pylith.testing.TestCases import TestAbstractComponent, make_suite
from pylith.problems.ProgressMonitor import ProgressMonitor


class TestProgressMonitor(TestAbstractComponent):
    """Unit testing of ProgressMonitor object.
    """
    _class = ProgressMonitor


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestProgressMonitor]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
