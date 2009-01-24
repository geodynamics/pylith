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
    io = MeshIOAscii()
    return


  def test_filename(self):
    """
    Test filename().
    """
    value = "hi.txt"

    io = MeshIOAscii()
    io.filename(value)
    self.assertEqual(value, io.filename())
    return


  def test_readwrite(self):
    """
    Test write() and read().
    """
    filenameIn = "data/mesh2Din3D.txt"
    filenameOut = "data/mesh2Din3D_test.txt"
    dim = 2

    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()
    cs._configure()

    io = MeshIOAscii()
    io.inventory.filename = filenameIn
    io.inventory.coordsys = cs
    io._configure()
    
    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()

    mesh = io.read(dim, normalizer, debug=False, interpolate=False)

    self.assertEqual(2, mesh.dimension())

    io.filename(filenameOut)
    io.write(mesh)

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
