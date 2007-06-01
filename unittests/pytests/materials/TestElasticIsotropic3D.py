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

## @file unittests/pytests/materials/TestElasticIsotropic3D.py

## @brief Unit testing of ElasticIsotropic3D object.

import unittest

from pylith.materials.ElasticIsotropic3D import ElasticIsotropic3D

# ----------------------------------------------------------------------
class TestElasticIsotropic3D(unittest.TestCase):
  """
  Unit testing of ElasticIsotropic3D object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    from pylith.materials.ElasticIsotropic3D import ElasticIsotropic3D
    material = ElasticIsotropic3D()

    self.assertNotEqual(None, material.cppHandle)
    return


  def test_dimension(self):
    """
    Test dimension().
    """
    material = ElasticIsotropic3D()
    self.assertEqual(3, material.dimension)
    return


  def test_useElasticBehavior(self):
    """
    Test useElasticBehavior().
    """
    material = ElasticIsotropic3D()
    material.useElasticBehavior(False)
    return


# End of file 
