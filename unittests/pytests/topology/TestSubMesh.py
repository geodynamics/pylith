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

## @file unittests/pytests/topology/TestSubMesh.py

## @brief Unit testing of Mesh object.

import unittest

from pylith.topology.Mesh import Mesh

# ----------------------------------------------------------------------
class TestSubMesh(unittest.TestCase):
  """
  Unit testing of Mesh object.
  """

  def test_constructorA(self):
    """
    Test constructor.
    """
    mesh = Mesh()
    self.assertEqual(0, mesh.dimension())
    self.assertEqual(False, mesh.debug())
    return


  def test_constructorB(self):
    """
    Test constructor.
    """
    mesh = self._getMesh()
    submesh = Mesh(mesh=mesh, label="bc")
    self.assertEqual(1, submesh.dimension())
    self.assertEqual(False, mesh.debug())
    return


  def test_coordsys(self):
    """
    Test coordsys().
    """
    mesh = self._getMesh()
    submesh = Mesh(mesh=mesh, label="bc")
    self.assertEqual(2, submesh.coordsys().spaceDim())
    return


  def test_debug(self):
    """
    Test debug().
    """
    mesh = self._getMesh()
    submesh = Mesh(mesh=mesh, label="bc")

    self.assertEqual(False, submesh.debug())

    submesh.debug(True)
    self.assertEqual(True, submesh.debug())
    return


  def test_dimension(self):
    """
    Test debug().
    """
    mesh = self._getMesh()
    submesh = Mesh(mesh=mesh, label="bc")

    self.assertEqual(1, submesh.dimension())
    return


  def test_comm(self):
    """
    Test comm().
    """
    mesh = self._getMesh()
    submesh = Mesh(mesh=mesh, label="bc")

    comm = submesh.comm()
    self.assertEqual(0, comm.rank)
    self.assertEqual(1, comm.size)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _getMesh(self):
    """
    Get mesh from file.
    """
    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()
    cs.inventory.spaceDim = 2
    cs._configure()

    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    normalizer._configure()    

    from pylith.meshio.MeshIOAscii import MeshIOAscii
    importer = MeshIOAscii()
    importer.inventory.filename = "data/tri3.mesh"
    importer.inventory.coordsys = cs
    importer._configure()
    mesh = importer.read(debug=False, interpolate=False)
    
    return mesh
  

# End of file 
