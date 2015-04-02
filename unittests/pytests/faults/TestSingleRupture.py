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

## @file unittests/pytests/faults/TestSingleRupture.py

## @brief Unit testing of SingleRupture object.

import unittest

# ----------------------------------------------------------------------
class TestSingleRupture(unittest.TestCase):
  """
  Unit testing of SingleRupture object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    from pylith.faults.SingleRupture import SingleRupture
    faults = SingleRupture()
    return


  def test_configure(self):
    """
    Test _configure().
    """
    from pylith.faults.SingleRupture import SingleRupture
    faults = SingleRupture()
    from pylith.faults.EqKinSrc import EqKinSrc
    faults.inventory.rupture = EqKinSrc()
    faults._configure()
    self.assertEqual(1, len(faults.components()))
    return


# End of file 
