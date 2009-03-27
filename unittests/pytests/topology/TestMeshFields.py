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

## @file unittests/pytests/topology/TestMeshFields.py

## @brief Unit testing of MeshFields object.

import unittest

from pylith.topology.Fields import MeshFields

# ----------------------------------------------------------------------
class TestMeshFields(unittest.TestCase):
  """
  Unit testing of MeshFields object.
  """

  def setUp(self):
    """
    Setup mesh and associated field.
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
    self.mesh = importer.read(normalizer, debug=False, interpolate=False)
    
    self.fields = MeshFields(self.mesh)
    return


  def test_constructorA(self):
    """
    Test constructor.
    """
    return


  def test_mesh(self):
    """
    Test mesh().
    """
    mesh = self.fields.mesh()
    
    self.assertEqual(2, mesh.dimension())
    return


  def test_add(self):
    self.fields.add("field", "displacement")
    field = self.fields.get("field")
    field.newSection()

    self.assertEqual(2, field.spaceDim())
    return


  def test_addFiberDim(self):
    from pylith.topology.topology import FieldBase
    self.fields.add("field", "displacement")
    field = self.fields.get("field")
    field.newSection(FieldBase.VERTICES_FIELD, 4)

    self.assertEqual(2, field.spaceDim())
    return


  def test_del(self):
    self.fields.add("field A", "A")
    self.fields.add("field B", "B")
    self.fields.delField("field A")
    field = self.fields.get("field B")
    return


  def test_copyLayout(self):
    from pylith.topology.topology import FieldBase
    self.fields.add("field A", "A")
    field = self.fields.get("field A")
    field.newSection(FieldBase.VERTICES_FIELD, 4)

    self.fields.add("field B", "B")
    self.fields.copyLayout("field A")

    # No test of result
    return


# End of file 
