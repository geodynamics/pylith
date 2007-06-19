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

## @file unittests/pytests/materials/TestElasticPlaneStrain.py

## @brief Unit testing of ElasticPlaneStrain object.

import unittest

from pylith.materials.ElasticPlaneStrain import ElasticPlaneStrain

# ----------------------------------------------------------------------
class TestElasticPlaneStrain(unittest.TestCase):
  """
  Unit testing of ElasticPlaneStrain object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    material = ElasticPlaneStrain()
    material._createCppHandle()
    self.assertNotEqual(None, material.cppHandle)
    return


  def test_dimension(self):
    """
    Test dimension().
    """
    material = ElasticPlaneStrain()
    material._createCppHandle()
    self.assertEqual(2, material.dimension)
    return


  def test_useElasticBehavior(self):
    """
    Test useElasticBehavior().
    """
    material = ElasticPlaneStrain()
    material._createCppHandle()
    material.useElasticBehavior(False)
    return


# End of file 
