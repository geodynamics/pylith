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

# ----------------------------------------------------------------------
class TestElasticPlaneStress(unittest.TestCase):
  """
  Unit testing of ElasticPlaneStress object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    from pylith.materials.ElasticPlaneStress import ElasticPlaneStress
    material = ElasticPlaneStress()

    self.assertNotEqual(None, material.cppHandle)
    return


# End of file 
