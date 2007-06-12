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

## @file unittests/pytests/bc/TestBCTwoSides.py

## @brief Unit testing of BCTwoSides object.

import unittest

# ----------------------------------------------------------------------
class TestBCTwoSides(unittest.TestCase):
  """
  Unit testing of BCTwoSides object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    from pylith.bc.BCTwoSides import BCTwoSides
    bc = BCTwoSides()
    return


  def test_configure(self):
    """
    Test _configure().
    """
    from pylith.bc.BCTwoSides import BCTwoSides
    bc = BCTwoSides()
    from pylith.bc.Dirichlet import Dirichlet
    bc.inventory.neg = Dirichlet()
    bc.inventory.pos = Dirichlet()
    bc._configure()
    self.assertEqual(2, len(bc.bin))
    return


# End of file 
