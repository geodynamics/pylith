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
