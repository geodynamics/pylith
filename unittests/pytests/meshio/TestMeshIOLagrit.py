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

## @file unittests/pytests/meshio/TestMeshIOLagrit.py

## @brief Unit testing of Python MeshIOLagrit object.

import unittest

from pylith.meshio.MeshIOLagrit import MeshIOLagrit
from pylith.meshio.MeshIOAscii import MeshIOAscii

# ----------------------------------------------------------------------
class TestMeshIOLagrit(unittest.TestCase):
  """
  Unit testing of Python MeshIOLagrit object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    iohandler = MeshIOLagrit()
    self.assertNotEqual(None, iohandler.cppHandle)
    return


  def test_filename(self):
    """
    Test filename().
    """
    iohandler = MeshIOLagrit()
    valueGmv = "hi.txt"
    valuePset = "hi2.txt"
    iohandler.filenameGmv = valueGmv
    iohandler.filenamePset = valuePset
    self.assertEqual(valueGmv, iohandler.filenameGmv)
    self.assertEqual(valuePset, iohandler.filenamePset)
    return


  def test_readwrite(self):
    """
    Test read().
    """
    from spatialdata.geocoords.CSCart import CSCart

    # For now, we only test reading the file. We would like to write
    # the file and compare against the original.
    iohandler = MeshIOLagrit()

    filenameGmvIn = "data/cube2_ascii.gmv"
    filenamePsetIn = "data/cube2_ascii.pset"
    filenameOut = "data/cube2_test.txt"
    filenameE = "data/cube2.txt"
    
    iohandler.filenameGmv = filenameGmvIn
    iohandler.filenamePset = filenamePsetIn
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
