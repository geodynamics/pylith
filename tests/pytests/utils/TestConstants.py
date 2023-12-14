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
# @file tests/pytests/utils/TestEventLogger.py

# @brief Unit testing of EventLogger object.

import unittest


# ----------------------------------------------------------------------
class TestConstants(unittest.TestCase):
  """Unit testing of constants.
  """
  

  def test_maxdouble(self):
    from pylith.utils.utils import maxdouble
    self.assertAlmostEqual(1.0, maxdouble()/1.0e+99, 7)
    return


  def test_maxfloat(self):
    from pylith.utils.utils import maxfloat
    self.assertAlmostEqual(1.0, maxfloat()/1.0e+30, 7)
    return


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestConstants))

    from pylith.utils.PetscManager import PetscManager
    petsc = PetscManager()
    petsc.initialize()

    success = unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful()

    petsc.finalize()


# End of file 
