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

## @file unittests/pytests/materials/TestElasticStress1D.py

## @brief Unit testing of ElasticStress1D object.

import unittest

# ----------------------------------------------------------------------
class TestElasticStress1D(unittest.TestCase):
  """
  Unit testing of ElasticStress1D object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    from pylith.materials.ElasticStress1D import ElasticStress1D
    material = ElasticStress1D()

    self.assertNotEqual(None, material.cppHandle)
    return


# End of file 
