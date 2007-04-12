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

# ----------------------------------------------------------------------
class TestElasticStrain1D(unittest.TestCase):
  """
  Unit testing of ElasticStrain1D object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    from pylith.materials.ElasticStrain1D import ElasticStrain1D
    material = ElasticStrain1D()

    self.assertNotEqual(None, material.cppHandle)
    return


# End of file 
