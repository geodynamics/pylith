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

## @file unittests/pytests/topology/TestMesh.py

## @brief Unit testing of Mesh object.

import unittest

from pylith.topology.Mesh import Mesh

# ----------------------------------------------------------------------
class TestMesh(unittest.TestCase):
  """
  Unit testing of Mesh object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    mesh = Mesh()
    self.assertNotEqual(None, mesh.cppHandle)
    return


  def test_initialize(self):
    """
    Test initialize().
    """
    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()

    mesh = Mesh()
    mesh.initialize(cs)
    self.assertEqual(cs, mesh.coordsys)
    return


# End of file 
