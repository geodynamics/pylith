#!/usr/bin/env nemesis
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#
# @file tests/pytests/utils/TestPetscManager.py
#
# @brief Unit testing of Python PetscManager object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.utils.PetscManager import (PetscManager, petsc_manager)


class TestPetscManager(TestComponent):
    """Unit testing of PetscManager object.
    """
    _class = PetscManager
    _factory = petsc_manager


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestPetscManager))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
