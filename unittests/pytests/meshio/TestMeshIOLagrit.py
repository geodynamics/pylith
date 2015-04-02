#!/usr/bin/env python
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
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
    io = MeshIOLagrit()
    return


  def test_filename(self):
    """
    Test filename().
    """
    valueGmv = "hi.txt"
    valuePset = "hi2.txt"

    io = MeshIOLagrit()
    io.filenameGmv(valueGmv)
    io.filenamePset(valuePset)
    self.assertEqual(valueGmv, io.filenameGmv())
    self.assertEqual(valuePset, io.filenamePset())
    return


  def test_readwrite(self):
    """
    Test read().
    """
    filenameGmvIn = "data/cube2_ascii.gmv"
    filenamePsetIn = "data/cube2_ascii.pset"
    filenameOut = "data/cube2_test.txt"
    filenameE = "data/cube2.txt"

    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()
    cs._configure()

    # For now, we only test reading the file. We would like to write
    # the file and compare against the original.
    io = MeshIOLagrit()
    io.inventory.filenameGmv = filenameGmvIn
    io.inventory.filenamePset = filenamePsetIn
    io._configure()

    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()

    mesh = io.read(debug=False, interpolate=True)

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


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.MeshIOLagrit import mesh_io
    io = mesh_io()
    return


# End of file 
