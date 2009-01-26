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

  def test_implementsConstraint(self):
    """
    Test to make sure DirichletBC satisfies constraint requirements.
    """
    bc = DirichletBC()
    from pylith.feassemble.Constraint import implementsConstraint
    self.failUnless(implementsConstraint(bc))
    return
    

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

    (mesh, bc, fields) = self._initialize()

    self.assertNotEqual(None, bc.cppHandle)

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return


  def test_setConstraintSizes(self):
    """
    Test setConstraintSizes().

    WARNING: This is not a rigorous test of setConstraintSizes() because we
    don't verify the results.
    """

    (mesh, bc, fields) = self._initialize()
    field = fields.getReal("field")
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

    (mesh, bc, fields) = self._initialize()
    field = fields.getReal("field")
    bc.setConstraintSizes(field)
    mesh.allocateRealSection(field)
    bc.setConstraints(field)

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return


  def test_useSolnIncr(self):
    """
    Test useSolnIncr().
    """
    (mesh, bc, fields) = self._initialize()
    bc.useSolnIncr(True)
    return


  def test_setField(self):
    """
    Test setField().

    WARNING: This is not a rigorous test of setField() because we
    don't verify the results.
    """

    (mesh, bc, fields) = self._initialize()
    field = fields.getReal("field")
    bc.setConstraintSizes(field)
    mesh.allocateRealSection(field)
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
    (mesh, bc, fields) = self._initialize()
    bc.finalize()

    # We should really add something here to check to make sure things
    # actually initialized correctly.
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _initialize(self):
    """
    Initialize DirichletBC boundary condition.
    """
    from pylith.bc.DirichletBC import DirichletBC
    bc = DirichletBC()
    bc._configure()
    bc.id = 0
    bc.label = "bc"
    bc.fixedDOF = [1]

    from pyre.units.time import second
    bc.tRef = -1.0*second

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    db = SimpleDB()
    db._configure()
    db.label = "TestDirichletBC tri3"
    db.iohandler.filename = "data/tri3.spatialdb"
    db.initialize()
    bc.db = db

    from pylith.bc.FixedDOFDB import FixedDOFDB
    dbRate = FixedDOFDB()
    dbRate._configure()
    dbRate.label = "TestDirichletBC rate tri3"
    dbRate.initialize()
    bc.dbRate = dbRate

    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()
    cs.spaceDim = 2

    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    normalizer.initialize()

    from pylith.meshio.MeshIOAscii import MeshIOAscii
    importer = MeshIOAscii()
    importer.filename = "data/tri3.mesh"
    importer.coordsys = cs
    mesh = importer.read(normalizer, debug=False, interpolate=False)
    
    bc.preinitialize(mesh)
    bc.initialize(totalTime=0.0, numTimeSteps=1, normalizer=normalizer)

    # Setup fields
    from pylith.topology.FieldsManager import FieldsManager
    fields = FieldsManager(mesh)
    fields.addReal("field")
    fields.setFiberDimension("field", cs.spaceDim)
    fields.allocate("field")

    import pylith.topology.topology as bindings
    bindings.zeroRealSection(fields.getReal("field"))
    
    return (mesh, bc, fields)


# End of file 
