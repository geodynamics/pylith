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

## @file unittests/pytests/materials/TestElasticStrain1D.py

## @brief Unit testing of ElasticStrain1D object.

import unittest

from pylith.materials.ElasticStrain1D import ElasticStrain1D

# ----------------------------------------------------------------------
class TestElasticStrain1D(unittest.TestCase):
  """
  Unit testing of ElasticStrain1D object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    material = ElasticStrain1D()

    self.assertEqual(1, material.dimension())
    return


  def test_useElasticBehavior(self):
    """
    Test useElasticBehavior().
    """
    material = ElasticStrain1D()
    material.useElasticBehavior(False)
    return


# End of file 
