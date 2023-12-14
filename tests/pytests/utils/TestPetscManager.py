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
# @file tests/pytests/utils/TestPetscManager.py
#
# @brief Unit testing of Python PetscManager object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.utils.PetscManager import (PetscManager, property_list)


class TestPetscManager(TestComponent):
    """Unit testing of PetscManager object.
    """
    _class = PetscManager
    _factory = property_list


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestPetscManager))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
