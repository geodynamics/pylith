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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
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


  def test_configure(self):
    """
    Test constructor.
    """
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    db = SimpleDB()
    db.inventory.label = "simple database"
    db._configure()

    from spatialdata.spatialdb.TimeHistory import TimeHistory
    th = TimeHistory()
    th._configure()

    from pylith.bc.DirichletBC import DirichletBC
    bc = DirichletBC()
    bc.inventory.label = "abc"
    bc.inventory.dbInitial = db
    bc.inventory.dbRate = db
    bc.inventory.dbChange = db
    bc.inventory.thChange = th    
    bc._configure()
    return


  def test_implementsConstraint(self):
    """
    Test to make sure DirichletBC satisfies constraint requirements.
    """
    bc = DirichletBC()
    from pylith.feassemble.Constraint import implementsConstraint
    self.failUnless(implementsConstraint(bc))
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


  def test_numDimConstrained(self):
    """
    Test numDimConstrained().
    """

    (mesh, bc, field) = self._initialize()

    self.assertEqual(1, bc.numDimConstrained())
    return


  def test_verifyConfiguration(self):
    """
    Test verifyConfiguration().

    WARNING: This is not a rigorous test of verifyConfiguration() because we
    don't verify the results.
    """

    (mesh, bc, field) = self._initialize()
    field.allocate()
    bc.verifyConfiguration()

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


  def test_setFieldIncr(self):
    """
    Test setFieldIncr().

    WARNING: This is not a rigorous test of setField() because we
    don't verify the results.
    """

    (mesh, bc, field) = self._initialize()
    bc.setConstraintSizes(field)
    field.allocate()
    bc.setConstraints(field)
    t0 = 1.0
    t1 = 2.0
    bc.setFieldIncr(t0, t1, field)

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


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.bc.DirichletBC import boundary_condition
    bc = boundary_condition()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _initialize(self):
    """
    Initialize DirichletBC boundary condition.
    """
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    dbInitial = SimpleDB()
    dbInitial.inventory.label = "TestDirichletBC tri3"
    dbInitial.inventory.iohandler.inventory.filename = "data/tri3_disp.spatialdb"
    dbInitial.inventory.iohandler._configure()
    dbInitial._configure()

    from pylith.bc.DirichletBC import DirichletBC
    bc = DirichletBC()
    bc.inventory.label = "bc"
    bc.inventory.bcDOF = [1]
    bc.inventory.dbInitial = dbInitial
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
    mesh = importer.read(debug=False, interpolate=False)
    
    bc.preinitialize(mesh)
    bc.initialize(totalTime=0.0, numTimeSteps=1, normalizer=normalizer)

    # Setup field
    from pylith.topology.Field import Field
    field = Field(mesh)
    field.newSection(field.VERTICES_FIELD, cs.spaceDim())
    
    return (mesh, bc, field)


# End of file 
