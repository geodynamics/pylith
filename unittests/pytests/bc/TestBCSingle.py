#!/usr/bin/env python
#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# {LicenseText}
#
# ======================================================================
#

## @file unittests/pytests/bc/TestBCSingle.py

## @brief Unit testing of BCSingle object.

import unittest

# ----------------------------------------------------------------------
class TestBCSingle(unittest.TestCase):
  """
  Unit testing of BCSingle object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    from pylith.bc.BCSingle import BCSingle
    bc = BCSingle()
    return


  def test_configure(self):
    """
    Test _configure().
    """
    from pylith.bc.BCSingle import BCSingle
    bc = BCSingle()
    from pylith.bc.Dirichlet import Dirichlet
    bc.inventory.bc = Dirichlet()
    bc._configure()
    self.assertEqual(1, len(bc.bin))
    return


# End of file 
