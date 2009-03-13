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

## @file unittests/pytests/materials/TestElasticPlaneStress.py

## @brief Unit testing of ElasticPlaneStress object.

import unittest

from pylith.materials.ElasticPlaneStress import ElasticPlaneStress

# ----------------------------------------------------------------------
class TestElasticPlaneStress(unittest.TestCase):
  """
  Unit testing of ElasticPlaneStress object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    material = ElasticPlaneStress()

    self.assertEqual(2, material.dimension())
    return


  def test_useElasticBehavior(self):
    """
    Test useElasticBehavior().
    """
    material = ElasticPlaneStress()
    material.useElasticBehavior(False)
    return


# End of file 
