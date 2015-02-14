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

## @file unittests/pytests/topology/TestMesh.py

## @brief Unit testing of Mesh object.

import unittest

from pylith.topology.Mesh import Mesh

# ----------------------------------------------------------------------
class TestMesh(unittest.TestCase):
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
    dim = 2
    mesh = Mesh(dim=dim)
    self.assertEqual(dim, mesh.dimension())
    self.assertEqual(False, mesh.debug())
    return


  def test_constructorC(self):
    """
    Test constructor.
    """
    dim = 2
    from pylith.mpi.Communicator import mpi_comm_self
    mesh = Mesh(dim=dim, comm=mpi_comm_self())
    self.assertEqual(dim, mesh.dimension())
    self.assertEqual(False, mesh.debug())
    return


  def test_coordsys(self):
    """
    Test coordsys().
    """
    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()
    cs.inventory.spaceDim = 2
    cs._configure()

    mesh = Mesh()
    mesh.coordsys(cs)
    self.assertEqual(cs.spaceDim(), mesh.coordsys().spaceDim())
    return


  def test_debug(self):
    """
    Test debug().
    """
    mesh = Mesh()

    self.assertEqual(False, mesh.debug())

    mesh.debug(True)
    self.assertEqual(True, mesh.debug())
    return


  def test_dimension(self):
    """
    Test debug().
    """
    mesh = Mesh()
    self.assertEqual(0, mesh.dimension())
    return


  def test_comm(self):
    """
    Test comm().
    """
    from pylith.mpi.Communicator import petsc_comm_self
    mesh = Mesh(dim=3, comm=petsc_comm_self())
    comm = mesh.comm()
    self.assertEqual(0, comm.rank)
    self.assertEqual(1, comm.size)
    return


  def test_view(self):
    """
    Test view().
    """
    mesh = self._getMesh()

    mesh.view()

    import os
    filename = "mesh.view"
    if os.path.isfile(filename):
      os.remove(filename)
    mesh.view(":%s:ascii_info_detail" % filename)
    self.assertTrue(os.path.isfile(filename))
    return


  def test_checkMaterialIds(self):
    """
    Test checkMaterialIds().
    """
    mesh = self._getMesh()
    mesh.checkMaterialIds([3, 4])
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
