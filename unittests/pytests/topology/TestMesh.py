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

  def test_constructor(self):
    """
    Test constructor.
    """
    mesh = Mesh()
    self.assertNotEqual(None, mesh.cppHandle)
    return


  def test_initialize(self):
    """
    Test initialize().
    """
    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()

    mesh = Mesh()
    mesh.initialize(cs)
    self.assertEqual(cs, mesh.coordsys)
    return


  def test_comm(self):
    """
    Test comm().
    """
    mesh = self._getMesh()
    comm = mesh.comm()
    import mpi
    self.assertEqual(mpi.MPI_COMM_WORLD, comm)
    return


  def test_getRealSection(self):
    """
    Test createRealSection().

    WARNING: This is not a rigorous test of createRealSection()
    because we don't verify the results.
    """
    mesh = self._getMesh()
    field = mesh.getRealSection("field")

    # We should really add something here to check to make sure things
    # actually initialized correctly.
    return


  def test_createRealSection(self):
    """
    Test createRealSection().

    WARNING: This is not a rigorous test of createRealSection()
    because we don't verify the results.
    """
    mesh = self._getMesh()
    field = mesh.createRealSection("field", fiberDim=2)

    # We should really add something here to check to make sure things
    # actually initialized correctly.
    return


  def test_allocateRealSection(self):
    """
    Test createRealSection().

    WARNING: This is not a rigorous test of createRealSection()
    because we don't verify the results.
    """
    mesh = self._getMesh()
    field = mesh.createRealSection("field", fiberDim=2)
    mesh.allocateRealSection(field)

    # We should really add something here to check to make sure things
    # actually initialized correctly.
    return


  def test_createMatrix(self):
    """
    Test createRealSection().

    WARNING: This is not a rigorous test of createRealSection()
    because we don't verify the results.
    """
    mesh = self._getMesh()
    field = mesh.createRealSection("field", fiberDim=2)
    mesh.allocateRealSection(field)
    matrix = mesh.createMatrix(field)

    # We should really add something here to check to make sure things
    # actually initialized correctly.
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _getMesh(self):
    """
    Get mesh from file.
    """
    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()
    cs.spaceDim = 2

    from pylith.meshio.MeshIOAscii import MeshIOAscii
    importer = MeshIOAscii()
    importer.filename = "data/tri3.mesh"
    importer.coordsys = cs
    mesh = importer.read(debug=False, interpolate=False)
    
    return mesh
  

# End of file 
