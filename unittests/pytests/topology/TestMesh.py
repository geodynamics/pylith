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


  def test_createSieveMesh(self):
    """
    Test createSeiveMesh().
    """
    mesh = Mesh()

    mesh.createSieveMesh()
    self.assertEqual(3, mesh.dimension())

    dim = 2
    mesh.createSieveMesh(dim)
    self.assertEqual(dim, mesh.dimension())
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
    mesh = Mesh()
    mesh.setComm(petsc_comm_self())
    comm = mesh.getComm()
    self.assertEqual(0, comm.rank)
    self.assertEqual(1, comm.size)
    return


  def test_view(self):
    """
    Test view().
    """
    mesh = self._getMesh()

    mesh.view("Testing view")
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
