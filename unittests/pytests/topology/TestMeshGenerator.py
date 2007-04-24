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

## @file unittests/pytests/topology/TestMeshGenerator.py

## @brief Unit testing of MeshGenerator object.

import unittest

from pylith.topology.MeshGenerator import MeshGenerator

# ----------------------------------------------------------------------
class TestMeshGenerator(unittest.TestCase):
  """
  Unit testing of MeshIO object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    generator = MeshGenerator()
    return
  

  def test_debug(self):
    """
    Test debug().
    """
    generator = MeshGenerator()

    value = False # default should be False
    self.assertTrue(value == generator.debug)

    value = True
    generator.debug = value
    self.assertTrue(value == generator.debug)
    return


  def test_interpolate(self):
    """
    Test interpolate access.
    """
    generator = MeshGenerator()

    value = False # default should be False
    self.assertEqual(value, generator.interpolate)

    value = True
    generator.interpolate = value
    self.assertEqual(value, generator.interpolate)
    return


  def test_distribution(self):
    """
    Test mesh distribution
    """
    generator = MeshGenerator()
    generator.interpolate = True
    generator.create(generator.createCubeBoundary())
    self.assertEqual(True, generator.interpolate)
    return


# End of file 
