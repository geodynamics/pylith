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

## @file unittests/pytests/materials/TestHomogeneous.py

## @brief Unit testing of Homogenous object.

import unittest

# ----------------------------------------------------------------------
class TestHomogeneous(unittest.TestCase):
  """
  Unit testing of Material object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    from pylith.materials.Homogeneous import Homogeneous
    materials = Homogeneous()
    return


# End of file 
