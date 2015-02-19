#!/usr/bin/env python
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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file unittests/pytests/utils/TestEventLogger.py

## @brief Unit testing of EventLogger object.

import unittest


# ----------------------------------------------------------------------
class TestConstants(unittest.TestCase):
  """
  Unit testing of constants.
  """
  

  def test_maxdouble(self):
    """
    Test maxdouble()
    """
    from pylith.utils.utils import maxdouble
    self.assertAlmostEqual(1.0, maxdouble()/1.0e+99, 7)
    return


  def test_maxflat(self):
    """
    Test maxflat()
    """
    from pylith.utils.utils import maxfloat
    self.assertAlmostEqual(1.0, maxfloat()/1.0e+30, 7)
    return


# End of file 
