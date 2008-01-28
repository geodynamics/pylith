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
    from pylith.bc.DirichletPoints import DirichletPoints
    bc.inventory.xNeg = DirichletPoints()
    bc.inventory.xPos = DirichletPoints()
    bc.inventory.yNeg = DirichletPoints()
    bc.inventory.yPos = DirichletPoints()
    bc._configure()
    self.assertEqual(4, len(bc.bin))
    return


# End of file 
