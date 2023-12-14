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
# @file tests/pytests/utils/TestProfiling.py
#
# @brief Unit testing of Python EmptyBin object.

import unittest

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


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestProfiling))

    from pylith.utils.PetscManager import PetscManager
    petsc = PetscManager()
    petsc.initialize()

    success = unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful()

    petsc.finalize()


# End of file
