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
# @file tests/pytests/apps/TestPetscApplication.py
#
# @brief Unit testing of Python PetscApplication object.

import unittest

from pylith.apps.PetscApplication import PetscApplication
from pylith.testing.UnitTestApp import configureComponent


class TestPetscApplication(unittest.TestCase):
    """Unit testing of PetscApplication object.
    """

    def test_constructor(self):
        app = PetscApplication()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestPetscApplication))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
