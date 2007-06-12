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

## @file unittests/pytests/materials/TestBiMaterial.py

## @brief Unit testing of BiMaterial object.

import unittest

# ----------------------------------------------------------------------
class TestBiMaterial(unittest.TestCase):
  """
  Unit testing of BiMaterial object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    from pylith.materials.BiMaterial import BiMaterial
    materials = BiMaterial()
    return


  def test_configure(self):
    """
    Test _configure().
    """
    from pylith.materials.BiMaterial import BiMaterial
    materials = BiMaterial()
    from pylith.materials.ElasticIsotropic3D import ElasticIsotropic3D
    materials.inventory.one = ElasticIsotropic3D()
    materials.inventory.two = ElasticIsotropic3D()
    materials._configure()
    self.assertEqual(2, len(materials.bin))
    return


# End of file 
