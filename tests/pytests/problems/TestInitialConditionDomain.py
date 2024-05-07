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

from pylith.testing.TestCases import TestComponent, make_suite
from pylith.problems.InitialConditionDomain import (InitialConditionDomain, initial_conditions)


class TestInitialConditionDomain(TestComponent):
    """Unit testing of InitialConditionDomain object.
    """
    _class = InitialConditionDomain
    _factory = initial_conditions


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestInitialConditionDomain]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
