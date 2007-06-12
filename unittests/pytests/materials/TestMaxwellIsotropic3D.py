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

## @file unittests/pytests/materials/TestMaxwellIsotropic3D.py

## @brief Unit testing of MaxwellIsotropic3D object.

import unittest

from pylith.materials.MaxwellIsotropic3D import MaxwellIsotropic3D

# ----------------------------------------------------------------------
class TestMaxwellIsotropic3D(unittest.TestCase):
  """
  Unit testing of MaxwellIsotropic3D object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    from pylith.materials.MaxwellIsotropic3D import MaxwellIsotropic3D
    material = MaxwellIsotropic3D()

    self.assertNotEqual(None, material.cppHandle)
    return


  def test_dimension(self):
    """
    Test dimension().
    """
    material = MaxwellIsotropic3D()
    self.assertEqual(3, material.dimension)
    return


  def test_useMaxwellBehavior(self):
    """
    Test useMaxwellBehavior().
    """
    material = MaxwellIsotropic3D()
    material.useElasticBehavior(False)
    return


# End of file 
