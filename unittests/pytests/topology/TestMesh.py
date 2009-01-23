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
    return


  def test_constructorB(self):
    """
    Test constructor.
    """
    from pylith.mpi.Communicator import mpi_comm_self
    mesh = Mesh(comm=mpi_comm_self())
    return


  def test_constructorC(self):
    """
    Test constructor.
    """
    from pylith.mpi.Communicator import mpi_comm_self
    mesh = Mesh(comm=mpi_comm_self(), dim=2)
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
    
    dim = 3
    self.assertEqual(dim, mesh.dimension())
    return


  def test_initialize(self):
    """
    Test initialize().
    """
    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()

    mesh = Mesh()
    mesh.coordsys(cs)
    mesh.initialize()
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
    mesh = importer.read(normalizer, debug=False, interpolate=False)
    
    return mesh
  

# End of file 
