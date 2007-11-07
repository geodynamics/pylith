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

## @file unittests/pytests/bc/TestAbsorbingDampers.py

## @brief Unit testing of AbsorbingDampers object.

import unittest

from pylith.bc.AbsorbingDampers import AbsorbingDampers

from pyre.units.time import second

# ----------------------------------------------------------------------
class TestAbsorbingDampers(unittest.TestCase):
  """
  Unit testing of AbsorbingDampers object.
  """

  def test_implementsIntegrator(self):
    """
    Test to make sure AbsorbingDampers satisfies constraint requirements.
    """
    bc = AbsorbingDampers()
    from pylith.feassemble.Integrator import implementsIntegrator
    self.failUnless(implementsIntegrator(bc))
    return
    

  def test_constructor(self):
    """
    Test constructor.
    """
    from pylith.bc.AbsorbingDampers import AbsorbingDampers
    bc = AbsorbingDampers()
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


  def test_timeStep(self):
    """
    Test timeStep().
    """
    (mesh, bc, fields) = self._initialize()

    dt = 0.25*second
    bc.timeStep(dt)
    return

  
  def test_stableTimeStep(self):
    """
    Test stableTimeStep().
    """
    (mesh, bc, fields) = self._initialize()

    dt = 0.24*second
    bc.timeStep(dt)
    self.assertEqual(dt, bc.stableTimeStep())
    return

  
  def test_needNewJacobian(self):
    """
    Test needNewJacobian().
    """
    (mesh, bc, fields) = self._initialize()
    self.assertEqual(True, bc.needNewJacobian())
    return

  
  def test_useSolnIncr(self):
    """
    Test useSolnIncr().
    """
    (mesh, bc, fields) = self._initialize()
    bc.useSolnIncr(True)
    return


  def test_integrateResidual(self):
    """
    Test integrateResidual().

    WARNING: This is not a rigorous test of integrateResidual() because we
    don't verify the results.
    """
    (mesh, bc, fields) = self._initialize()

    residual = fields.getReal("residual")
    t = 0.02*second
    bc.integrateResidual(residual, t, fields)

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return


  def test_integrateJacobian(self):
    """
    Test integrateJacobian().

    WARNING: This is not a rigorous test of integrateJacobian() because we
    don't verify the results.
    """

    (mesh, bc, fields) = self._initialize()

    jacobian = mesh.createMatrix(fields.getReal("residual"))
    import pylith.utils.petsc as petsc
    petsc.mat_setzero(jacobian)
    t = 0.24*second
    bc.integrateJacobian(jacobian, t, fields)
    self.assertEqual(False, bc.needNewJacobian())

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return


  def test_updateState(self):
    """
    Test updateState().

    WARNING: This is not a rigorous test of updateState() because we
    neither set the input fields or verify the results.
    """
    (mesh, bc, fields) = self._initialize()

    dispT = fields.getReal("dispT")
    t = 0.50*second
    bc.updateState(t, dispT)

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
    Initialize AbsorbingDampers boundary condition.
    """
    from pylith.bc.AbsorbingDampers import AbsorbingDampers
    bc = AbsorbingDampers()
    bc._configure()
    bc.id = 0
    bc.label = "bc"

    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATSimplex()
    cell.shape = "line"
    cell.degree = 1
    cell.order = 1
    from pylith.feassemble.quadrature.Quadrature1Din2D import Quadrature1Din2D
    quadrature = Quadrature1Din2D()
    quadrature.cell = cell
    bc.quadrature = quadrature

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    db = SimpleDB()
    db._configure()
    db.label = "TestAbsorbingDampers tri3"
    db.iohandler.filename = "data/elasticplanestrain.spatialdb"
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
    
    bc.preinitialize(mesh)
    bc.initialize()
    bc.timeStep(0.01*second)

    # Setup fields
    from pylith.topology.FieldsManager import FieldsManager
    fields = FieldsManager(mesh)
    fields.addReal("residual")
    fields.addReal("solution")
    fields.addReal("dispT")
    fields.addReal("dispTmdt")
    fields.createHistory(["solution", "dispT", "dispTmdt"])
    fields.setFiberDimension("residual", cs.spaceDim)
    fields.allocate("residual")
    fields.copyLayout("residual")

    import pylith.topology.topology as bindings
    bindings.zeroRealSection(fields.getReal("residual"))
    
    return (mesh, bc, fields)


# End of file 
