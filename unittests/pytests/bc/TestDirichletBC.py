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

## @file unittests/pytests/bc/TestDirichletBC.py

## @brief Unit testing of DirichletBC object.

import unittest

from pylith.bc.DirichletBC import DirichletBC

# ----------------------------------------------------------------------
class TestDirichletBC(unittest.TestCase):
  """
  Unit testing of DirichletBC object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    from pylith.bc.DirichletBC import DirichletBC
    bc = DirichletBC()
    return


  def test_initialize(self):
    """
    Test initialize().

    WARNING: This is not a rigorous test of initialize() because we
    don't verify the results.
    """

    (mesh, bc, field) = self._initialize()

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return


  def test_setConstraintSizes(self):
    """
    Test setConstraintSizes().

    WARNING: This is not a rigorous test of setConstraintSizes() because we
    don't verify the results.
    """

    (mesh, bc, field) = self._initialize()

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

    (mesh, bc, field) = self._initialize()
    bc.setConstraintSizes(field)
    field.allocate()
    bc.setConstraints(field)

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return


  def test_useSolnIncr(self):
    """
    Test useSolnIncr().
    """
    (mesh, bc, field) = self._initialize()
    bc.useSolnIncr(True)
    return


  def test_setField(self):
    """
    Test setField().

    WARNING: This is not a rigorous test of setField() because we
    don't verify the results.
    """

    (mesh, bc, field) = self._initialize()
    bc.setConstraintSizes(field)
    field.allocate()
    bc.setConstraints(field)
    t = 1.0
    bc.setField(t, field)

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return


  def test_finalize(self):
    """
    Test finalize().

    WARNING: This is not a rigorous test of finalize() because we
    neither set the input fields or verify the results.
    """
    (mesh, bc, field) = self._initialize()
    bc.finalize()

    # We should really add something here to check to make sure things
    # actually initialized correctly.
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _initialize(self):
    """
    Initialize DirichletBC boundary condition.
    """
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    db = SimpleDB()
    db.inventory.label = "TestDirichletBC tri3"
    db.inventory.iohandler.inventory.filename = "data/tri3.spatialdb"
    db.inventory.iohandler._configure()
    db._configure()

    from pylith.bc.FixedDOFDB import FixedDOFDB
    dbRate = FixedDOFDB()
    dbRate.inventory.label = "TestDirichletBC rate tri3"
    dbRate._configure()

    from pylith.bc.DirichletBC import DirichletBC
    bc = DirichletBC()
    bc.inventory.label = "bc"
    bc.inventory.fixedDOF = [1]
    from pyre.units.time import second
    bc.inventory.tRef = -1.0*second
    bc.inventory.db = db
    bc.inventory.dbRate = dbRate
    bc._configure()

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
    
    bc.preinitialize(mesh)
    bc.initialize(totalTime=0.0, numTimeSteps=1, normalizer=normalizer)

    # Setup field
    from pylith.topology.Field import MeshField
    field = MeshField(mesh)
    field.newSection(field.VERTICES_FIELD, cs.spaceDim())
    field.allocate()

    field.zero()
    
    return (mesh, bc, field)


# End of file 
