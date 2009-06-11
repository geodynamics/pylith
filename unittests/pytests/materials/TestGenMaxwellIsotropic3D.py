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

## @file unittests/pytests/materials/TestGenMaxwellIsotropic3D.py

## @brief Unit testing of GenMaxwellIsotropic3D object.

import unittest

from pylith.materials.GenMaxwellIsotropic3D import GenMaxwellIsotropic3D

# ----------------------------------------------------------------------
class TestGenMaxwellIsotropic3D(unittest.TestCase):
  """
  Unit testing of GenMaxwellIsotropic3D object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    from pylith.materials.GenMaxwellIsotropic3D import GenMaxwellIsotropic3D
    material = GenMaxwellIsotropic3D()
    material._createCppHandle()
    self.assertNotEqual(None, material.cppHandle)
    return


  def test_dimension(self):
    """
    Test dimension().
    """
    material = GenMaxwellIsotropic3D()
    material._createCppHandle()
    self.assertEqual(3, material.dimension)
    return


  def test_useElasticBehavior(self):
    """
    Test useElasticBehavior().
    """
    material = GenMaxwellIsotropic3D()
    material._createCppHandle()
    material.useElasticBehavior(False)
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.materials.GenMaxwellIsotropic3D import material
    m = material()
    return


# End of file 
