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

## @file unittests/pytests/topology/TestMeshImporter.py

## @brief Unit testing of MeshImporter object.

import unittest

from pylith.topology.MeshImporter import MeshImporter

# ----------------------------------------------------------------------
class TestMeshImporter(unittest.TestCase):
  """
  Unit testing of MeshIO object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    importer = MeshImporter()
    return


# End of file 
