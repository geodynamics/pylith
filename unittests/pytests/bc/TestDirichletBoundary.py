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

## @file unittests/pytests/bc/TestDirichletBoundary.py

## @brief Unit testing of DirichletBoundary object.

import unittest

from pylith.bc.DirichletBoundary import DirichletBoundary

# ----------------------------------------------------------------------
class TestDirichletBoundary(unittest.TestCase):
  """
  Unit testing of DirichletBoundary object.
  """

  def test_implementsConstraint(self):
    """
    Test to make sure DirichletBoundary satisfies constraint requirements.
    """
    bc = DirichletBoundary()
    from pylith.feassemble.Constraint import implementsConstraint
    self.failUnless(implementsConstraint(bc))
    return
    

  def test_constructor(self):
    """
    Test constructor.
    """
    from pylith.bc.DirichletBoundary import DirichletBoundary
    bc = DirichletBoundary()
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
    from pyre.units.time import second
    t = 1.0*second
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
    Initialize DirichletBoundary boundary condition.
    """
    from pylith.bc.DirichletBoundary import DirichletBoundary
    bc = DirichletBoundary()
    bc._configure()
    bc.output._configure()
    bc.output.writer._configure()
    bc.label = "bc"
    bc.fixedDOF = [1]

    from pyre.units.time import second
    bc.tRef = -1.0*second

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    db = SimpleDB()
    db._configure()
    db.label = "TestDirichletBoundary tri3"
    db.iohandler.filename = "data/tri3.spatialdb"
    db.initialize()
    bc.db = db

    from pylith.bc.FixedDOFDB import FixedDOFDB
    dbRate = FixedDOFDB()
    dbRate._configure()
    dbRate.label = "TestDirichletBoundary rate tri3"
    dbRate.initialize()
    bc.dbRate = dbRate

    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()
    cs.spaceDim = 2

    from pylith.meshio.MeshIOAscii import MeshIOAscii
    importer = MeshIOAscii()
    importer.filename = "data/tri3.mesh"
    importer.coordsys = cs
    mesh = importer.read(debug=False, interpolate=False)
    
    bc.preinitialize(mesh)
    from pyre.units.time import second
    bc.initialize(totalTime=0.0*second, numTimeSteps=1)

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
