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

# ----------------------------------------------------------------------
class TestElasticPlaneStrain(unittest.TestCase):
  """
  Unit testing of ElasticPlaneStrain object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    from pylith.materials.ElasticPlaneStrain import ElasticPlaneStrain
    material = ElasticPlaneStrain()

    self.assertNotEqual(None, material.cppHandle)
    return


# End of file 
