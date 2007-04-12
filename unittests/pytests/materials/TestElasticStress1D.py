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

from pylith.materials.ElasticStress1D import ElasticStress1D

# ----------------------------------------------------------------------
class TestElasticStress1D(unittest.TestCase):
  """
  Unit testing of ElasticStress1D object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    material = ElasticStress1D()

    self.assertNotEqual(None, material.cppHandle)
    return


  def test_dimension(self):
    """
    Test dimension().
    """
    material = ElasticStress1D()
    self.assertEqual(1, material.dimension)
    return


# End of file 
