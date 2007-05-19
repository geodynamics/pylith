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

## @file unittests/pytests/topology/TestFieldsManager.py

## @brief Unit testing of FieldsManager object.

import unittest

from pylith.topology.FieldsManager import FieldsManager

# ----------------------------------------------------------------------
class TestFieldsManager(unittest.TestCase):
  """
  Unit testing of FieldsManager object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    mesh = self._initialize()
    from pylith.topology.FieldsManager import FieldsManager
    manager = FieldsManager(mesh)

    self.assertNotEqual(None, manager.cppHandle)
    return


  def test_addReal(self):
    """
    Test addReal().

    WARNING: This is not a rigorous test of initialize() because we
    don't verify the results.
    """

    mesh = self._initialize()
    from pylith.topology.FieldsManager import FieldsManager
    manager = FieldsManager(mesh)
    manager.addReal("field")

    # We should really add something here to check to make sure things
    # actually initialized correctly.
    return


  def test_getReal(self):
    """
    Test getReal().

    WARNING: This is not a rigorous test of setConstraintSizes() because we
    don't verify the results.
    """
    mesh = self._initialize()
    from pylith.topology.FieldsManager import FieldsManager
    manager = FieldsManager(mesh)
    manager.addReal("field")

    field = manager.getReal("field")

    # We should really add something here to check to make sure things
    # actually initialized correctly.
    return


  def test_delReal(self):
    """
    Test delReal().

    WARNING: This is not a rigorous test of setConstraints() because we
    don't verify the results.
    """

    mesh = self._initialize()
    from pylith.topology.FieldsManager import FieldsManager
    manager = FieldsManager(mesh)
    manager.addReal("field")

    manager.delReal("field")

    # We should really add something here to check to make sure things
    # actually initialized correctly.
    return


  def test_copyLayout(self):
    """
    Test copyLayout().

    WARNING: This is not a rigorous test of setField() because we
    don't verify the results.
    """

    mesh = self._initialize()
    from pylith.topology.FieldsManager import FieldsManager
    manager = FieldsManager(mesh)
    manager.addReal("field")

    # We should really add something here to check to make sure things
    # actually initialized correctly.
    return


  def test_copyLayoutFromSrc(self):
    """
    Test copyLayoutFromSrc().

    WARNING: This is not a rigorous test of setField() because we
    don't verify the results.
    """

    mesh = self._initialize()
    field = mesh.createRealSection("field", fiberDim=3)
    mesh.allocateRealSection(field)

    from pylith.topology.FieldsManager import FieldsManager
    manager = FieldsManager(mesh)
    manager.addReal("fieldA")
    manager.addReal("fieldB")

    manager.copyLayoutFromSrc(field)
    
    # We should really add something here to check to make sure things
    # actually initialized correctly.
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _initialize(self):
    """
    Initialize mesh.
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
