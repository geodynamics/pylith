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

## @file unittests/pytests/bc/TestBCFourSides.py

## @brief Unit testing of BCFourSides object.

import unittest

# ----------------------------------------------------------------------
class TestBCFourSides(unittest.TestCase):
  """
  Unit testing of BCFourSides object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    from pylith.bc.BCFourSides import BCFourSides
    bc = BCFourSides()
    return


  def test_configure(self):
    """
    Test _configure().
    """
    from pylith.bc.BCFourSides import BCFourSides
    bc = BCFourSides()
    from pylith.bc.Dirichlet import Dirichlet
    bc.inventory.xNeg = Dirichlet()
    bc.inventory.xPos = Dirichlet()
    bc.inventory.yNeg = Dirichlet()
    bc.inventory.yPos = Dirichlet()
    bc._configure()
    self.assertEqual(4, len(bc.bin))
    return


# End of file 
