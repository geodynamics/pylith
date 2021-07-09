#!/usr/bin/env nemesis
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ======================================================================
#
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
