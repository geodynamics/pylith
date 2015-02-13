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

    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()
    cs._configure()

    io = MeshIOAscii()
    io.inventory.filename = filenameIn
    io.inventory.coordsys = cs
    io._configure()
    
    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()

    mesh = io.read(debug=False, interpolate=True)

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


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.MeshIOAscii import mesh_io
    io = mesh_io()
    return


# End of file 
