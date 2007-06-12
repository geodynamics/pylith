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

## @file unittests/pytests/bc/TestBCSixSides.py

## @brief Unit testing of BCSixSides object.

import unittest

# ----------------------------------------------------------------------
class TestBCSixSides(unittest.TestCase):
  """
  Unit testing of BCSixSides object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    from pylith.bc.BCSixSides import BCSixSides
    bc = BCSixSides()
    return


  def test_configure(self):
    """
    Test _configure().
    """
    from pylith.bc.BCSixSides import BCSixSides
    bc = BCSixSides()
    from pylith.bc.Dirichlet import Dirichlet
    bc.inventory.xNeg = Dirichlet()
    bc.inventory.xPos = Dirichlet()
    bc.inventory.yNeg = Dirichlet()
    bc.inventory.yPos = Dirichlet()
    bc.inventory.zNeg = Dirichlet()
    bc.inventory.zPos = Dirichlet()
    bc._configure()
    self.assertEqual(6, len(bc.bin))
    return


# End of file 
