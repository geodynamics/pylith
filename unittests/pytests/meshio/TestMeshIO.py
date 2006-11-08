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

## @file unittests/pytests/meshio/TestMeshIO.py

## @brief Unit testing of MeshIO object.

import unittest

# ----------------------------------------------------------------------
class TestMeshIO(unittest.TestCase):
  """
  Unit testing of MeshIO object.
  """
  

  def test_interpolate(self):
    """
    Test interpolate access.
    """
    from pylith.meshio.MeshIO import MeshIO
    iohandler = MeshIO()

    value = False # default is False
    self.assertEqual(value, iohandler.interpolate)

    value = True
    iohandler.interpolate = value
    self.assertEqual(value, iohandler.interpolate)
    return


# End of file 
