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

## @file unittests/pytests/bc/TestDirichlet.py

## @brief Unit testing of Dirichlet object.

import unittest

from pylith.bc.Dirichlet import Dirichlet

# ----------------------------------------------------------------------
class TestDirichlet(unittest.TestCase):
  """
  Unit testing of Dirichlet object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    from pylith.bc.Dirichlet import Dirichlet
    bc = Dirichlet()

    self.assertNotEqual(None, bc.cppHandle)
    return


  def test_initialize(self):
    """
    Test initialize().

    WARNING: This is not a rigorous test of initialize() because we
    don't verify the results.
    """

    self._initialize()

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return


  def test_setConstraintSizes(self):
    """
    Test setConstraintSizes().

    WARNING: This is not a rigorous test of setConstraintSizes() because we
    don't verify the results.
    """

    (mesh, bc) = self._initialize()
    field = mesh.createRealSection("field", mesh.dimension())
    bc.setConstraintSizes(field)

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return


  def test_setConstraints(self):
    """
    Test setConstraints().

    WARNING: This is not a rigorous test of setConstraints() because we
    don't verify the results.
    """

    (mesh, bc) = self._initialize()
    field = mesh.createRealSection("field", mesh.dimension())
    bc.setConstraintSizes(field)
    mesh.allocateRealSection(field)
    bc.setConstraints(field)

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return


  def test_setField(self):
    """
    Test setField().

    WARNING: This is not a rigorous test of setField() because we
    don't verify the results.
    """

    (mesh, bc) = self._initialize()
    field = mesh.createRealSection("field", mesh.dimension())
    bc.setConstraintSizes(field)
    mesh.allocateRealSection(field)
    bc.setConstraints(field)
    t = 1.0
    bc.setField(t, field)

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _initialize(self):
    """
    Initialize Dirichlet boundary condition.
    """
    from pylith.bc.Dirichlet import Dirichlet
    bc = Dirichlet()
    bc.id = 0
    bc.label = "bc"
    bc.fixedDOF = [1]

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    iohandler = SimpleIOAscii()
    iohandler.filename = "data/tri3.spatialdb"
    db = SimpleDB()
    db.label = "TestDirichlet tri3"
    db.iohandler = iohandler
    db.initialize()
    bc.db = db

    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()
    cs.spaceDim = 2

    from pylith.meshio.MeshIOAscii import MeshIOAscii
    importer = MeshIOAscii()
    importer.filename = "data/tri3.mesh"
    importer.coordsys = cs
    mesh = importer.read(debug=False, interpolate=False)
    
    bc.initialize(mesh)
    return (mesh, bc)


# End of file 
