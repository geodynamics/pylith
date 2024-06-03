# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

import unittest

from pylith.testing.TestCases import make_suite
from pylith.utils.profiling import (resourceUsage, resourceUsageString)


class TestProfiling(unittest.TestCase):
    """Unit testing of profiling module.
    """

    def test_resourceUsage(self):
        cputime, memory = resourceUsage()
        self.assertFalse(cputime is None)
        self.assertFalse(memory is None)

    def test_resourceUsageString(self):
        usage = resourceUsageString()
        self.assertFalse(usage is None)


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestProfiling]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    from pylith.utils.PetscManager import PetscManager
    petsc = PetscManager()
    petsc.initialize()

    unittest.main(verbosity=2)

    petsc.finalize()


# End of file
