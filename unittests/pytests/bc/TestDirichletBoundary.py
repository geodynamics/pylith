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
    from pylith.feassemble.Constraint import Constraint
    self.failUnless(isinstance(bc, Constraint))
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
    Initialize DirichletBoundary boundary condition.
    """
    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()
    cs.spaceDim = 2

    from pylith.meshio.MeshIOAscii import MeshIOAscii
    importer = MeshIOAscii()
    importer.inventory.filename = "data/tri3.mesh"
    importer.inventory.coordsys = cs
    importer._configure()
    mesh = importer.read(normalizer, debug=False, interpolate=False)
    
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    db = SimpleDB()
    db.inventory.label = "TestDirichletBoundary tri3"
    db.inventory.iohandler.inventory.filename = "data/tri3.spatialdb"
    db.inventory.iohandler._configure()
    db._configure()

    from pylith.bc.FixedDOFDB import FixedDOFDB
    dbRate = FixedDOFDB()
    dbRate.inventory.label = "TestDirichletBoundary rate tri3"
    dbRate._configure()

    bc.inventory.db = db
    bc.inventory.dbRate = dbRate


    from pylith.bc.DirichletBoundary import DirichletBoundary
    bc = DirichletBoundary()
    bc.inventory.output._configure()
    bc.output.writer._configure()
    bc.label = "bc"
    bc.fixedDOF = [1]
    bc._configure()

    from pyre.units.time import second
    bc.tRef = -1.0*second

    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    normalizer.initialize()

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
    
    return (mesh, bc, field)


# End of file 
