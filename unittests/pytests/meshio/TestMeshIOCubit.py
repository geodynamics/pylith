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

## @file unittests/pytests/meshio/TestMeshIOCubit.py

## @brief Unit testing of Python MeshIOCubit object.

import unittest

from pylith.meshio.MeshIOCubit import MeshIOCubit
from pylith.meshio.MeshIOAscii import MeshIOAscii

# ----------------------------------------------------------------------
class TestMeshIOCubit(unittest.TestCase):
  """
  Unit testing of Python MeshIOCubit object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    iohandler = MeshIOCubit()
    self.assertNotEqual(None, iohandler.cppHandle)
    return


  def test_filename(self):
    """
    Test filename().
    """
    iohandler = MeshIOCubit()
    value = "hi.txt"
    iohandler.filename = value
    self.assertEqual(value, iohandler.filename)
    return


  def test_readwrite(self):
    """
    Test read().
    """
    from spatialdata.geocoords.CSCart import CSCart

    # For now, we only test reading the file.
    iohandler = MeshIOCubit()

    filenameIn = "data/twohex8.exo"
    filenameOut = "data/twohex8_test.txt"
    filenameE = "data/twohex8.txt"
    
    iohandler.filename = filenameIn
    iohandler.coordsys = CSCart()
    mesh = iohandler.read(debug=False, interpolate=False)
    self.assertEqual(3, mesh.dimension())

    testhandler = MeshIOAscii()
    testhandler.filename = filenameOut
    testhandler.coordsys = CSCart()
    testhandler.write(mesh)

    fileE = open(filenameE, "r")
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
