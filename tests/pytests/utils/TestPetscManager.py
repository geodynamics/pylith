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
from pylith.utils.PetscManager import (PetscManager, property_list)


class TestPetscManager(TestComponent):
    """Unit testing of PetscManager object.
    """
    _class = PetscManager
    _factory = property_list


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestPetscManager]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
