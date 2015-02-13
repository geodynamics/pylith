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
    field.allocate()
    bc.finalize()

    # We should really add something here to check to make sure things
    # actually initialized correctly.
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.bc.DirichletBoundary import boundary_condition
    bc = boundary_condition()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _initialize(self):
    """
    Initialize DirichletBoundary boundary condition.
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
    mesh = importer.read(debug=False, interpolate=False)
    
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    db = SimpleDB()
    db.inventory.label = "TestDirichletBoundary tri3"
    db.inventory.iohandler.inventory.filename = "data/tri3_disp.spatialdb"
    db.inventory.iohandler._configure()
    db._configure()

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    dbRate = SimpleDB()
    dbRate.inventory.label = "TestDirichletBoundary tri3"
    dbRate.inventory.iohandler.inventory.filename = "data/tri3_vel.spatialdb"
    dbRate.inventory.iohandler._configure()
    dbRate._configure()

    from pylith.bc.DirichletBoundary import DirichletBoundary
    bc = DirichletBoundary()
    bc.inventory.output._configure()
    bc.inventory.output.writer._configure()
    bc.inventory.label = "bc"
    bc.inventory.bcDOF = [1]
    bc.inventory.dbInitial = db
    bc.inventory.dbRate = dbRate
    bc._configure()

    bc.preinitialize(mesh)
    bc.initialize(totalTime=0.0, numTimeSteps=1, normalizer=normalizer)

    # Setup field
    from pylith.topology.Field import Field
    field = Field(mesh)
    field.newSection(field.VERTICES_FIELD, cs.spaceDim())
    
    return (mesh, bc, field)


# End of file 
