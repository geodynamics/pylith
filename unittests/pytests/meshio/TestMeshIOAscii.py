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


  def test_readwrite(self):
    """
    Test write() and read().
    """
    filenameIn = "data/mesh2Din3D.txt"
    filenameOut = "data/mesh2Din3D_test.txt"
    
    iohandler = MeshIOAscii()
    iohandler.filename = filenameIn
    mesh = iohandler.read()
    iohandler.filename = filenameOut
    iohandler.write(mesh)

    fileE = open(filenameIn, "r")
    linesE = fileE.readlines()
    fileE.close()
    fileT = open(filenameOut, "r")
    linesT = fileT.readlines()
    fileT.close()

    self.assertEqual(len(linesE), len(linesT))
    for (lineE, lineT) in zip(linesE, linesT):
      self.assertEqual(lineE, lineT)
    return


# End of file 
