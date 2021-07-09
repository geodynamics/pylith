#!/usr/bin/env python
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

## @file tests/pytests/utils/TestEventLogger.py

## @brief Unit testing of EventLogger object.

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
