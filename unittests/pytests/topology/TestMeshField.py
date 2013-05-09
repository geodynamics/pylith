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
# Copyright (c) 2010-2013 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file unittests/pytests/topology/TestMeshField.py

## @brief Unit testing of MeshField object.

import unittest

from pylith.topology.Field import MeshField

# ----------------------------------------------------------------------
class TestMeshField(unittest.TestCase):
  """
  Unit testing of MeshField object.
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
    self.mesh = importer.read(debug=False, interpolate=False)
    
    self.field = MeshField(self.mesh)
    self.field.allocate()
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
    mesh = self.field.mesh()
    
    self.assertEqual(2, mesh.dimension())
    return


  def test_label(self):
    """
    Test label().
    """
    label = "field A"

    self.field.label(label)
    self.assertEqual(label, self.field.label())
    return


  def test_vectorFieldType(self):
    """
    Test vectorFieldType().
    """
    fieldType = MeshField.MULTI_SCALAR

    self.field.vectorFieldType(fieldType)
    self.assertEqual(fieldType, self.field.vectorFieldType())
    return


  def test_scale(self):
    """
    Test scale().
    """
    scale = 2.0

    self.field.scale(scale)
    self.assertEqual(scale, self.field.scale())
    return


  def test_addDimensionOkay(self):
    """
    Test addDimensionOkay().
    """
    self.assertEqual(False, self.field.addDimensionOkay())

    self.field.addDimensionOkay(True)
    self.assertEqual(True, self.field.addDimensionOkay())
    return


  def test_spaceDim(self):
    """
    Test spaceDim().
    """
    self.assertEqual(2, self.field.spaceDim())
    return


  def test_newSection(self):
    """
    Test newSection().
    """
    self.field.newSection()

    # No test of result
    return


  def test_newSectionDomain(self):
    """
    Test newSection(domain).
    """
    self.field.newSection(MeshField.VERTICES_FIELD, 4)

    # No test of result
    return


  def test_cloneSectionField(self):
    """
    Test newSection(field).
    """
    fieldB = MeshField(self.mesh)
    fieldB.cloneSection(self.field)

    # No test of result
    return

  def test_operatorAdd(self):
    """
    Test add().
    """
    fieldB = MeshField(self.mesh)
    fieldB.allocate()
    self.field.add(fieldB)

    # No test of result
    return


  def test_copy(self):
    """
    Test newSection(field).
    """
    fieldB = MeshField(self.mesh)
    fieldB.allocate()
    fieldB.copy(self.field)

    # No test of result
    return


# End of file 
