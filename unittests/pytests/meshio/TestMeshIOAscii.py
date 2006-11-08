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

## @file unittests/pytests/meshio/TestMeshIOAscii.py

## @brief Unit testing of Python MeshIOAscii object.

import unittest
from pylith.meshio.MeshIOAscii import MeshIOAscii

# ----------------------------------------------------------------------
class TestMeshIOAscii(unittest.TestCase):
  """
  Unit testing of Python MeshIOAscii object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    iohandler = MeshIOAscii()
    self.failIf(None == iohandler.cppHandle)
    return


  def test_filename(self):
    """
    Test filename().
    """
    filename = "hi.txt"
    iohandler = MeshIOAscii()
    iohandler.filename = filename
    self.assertEqual(filename, iohandler.filename)
    return


  def test_writeread(self):
    """
    Test write() and read().
    """

    self.assertEqual(0,1)
    return


# End of file 
