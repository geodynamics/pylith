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
    from pylith.bc.DirichletPoints import DirichletPoints
    bc.inventory.xNeg = DirichletPoints()
    bc.inventory.xPos = DirichletPoints()
    bc.inventory.yNeg = DirichletPoints()
    bc.inventory.yPos = DirichletPoints()
    bc.inventory.zNeg = DirichletPoints()
    bc.inventory.zPos = DirichletPoints()
    bc._configure()
    self.assertEqual(6, len(bc.bin))
    return


# End of file 
