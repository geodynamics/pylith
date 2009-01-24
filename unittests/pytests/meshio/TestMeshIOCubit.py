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
    io = MeshIOCubit()
    return


  def test_filename(self):
    """
    Test filename().
    """
    value = "hi.txt"

    io = MeshIOCubit()
    io.filename(value)
    self.assertEqual(value, io.filename())
    return


  def test_readwrite(self):
    """
    Test read().
    """
    filenameIn = "data/twohex8.exo"
    filenameOut = "data/twohex8_test.txt"
    filenameE = "data/twohex8.txt"
    dim = 3

    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()
    cs._configure()

    # For now, we only test reading the file.
    io = MeshIOCubit()
    io.inventory.filename = filenameIn
    io._configure()
    
    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()

    mesh = io.read(dim, normalizer, debug=False, interpolate=False)
    self.assertEqual(3, mesh.dimension())

    testhandler = MeshIOAscii()
    testhandler.filename(filenameOut)
    testhandler.coordsys = cs
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
