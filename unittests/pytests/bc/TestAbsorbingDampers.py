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

# ----------------------------------------------------------------------
class TestAbsorbingDampers(unittest.TestCase):
  """
  Unit testing of AbsorbingDampers object.
  """

  def test_implementsIntegrator(self):
    """
    Test to make sure AbsorbingDampers satisfies integrator requirements.
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

    # No testing of result.
    return


  def test_timeStep(self):
    """
    Test timeStep().
    """
    (mesh, bc, fields) = self._initialize()

    dt = 0.25
    bc.timeStep(dt)
    return

  
  def test_stableTimeStep(self):
    """
    Test stableTimeStep().
    """
    (mesh, bc, fields) = self._initialize()

    self.assertEqual(1.0e+30, bc.stableTimeStep(mesh))
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

    residual = fields.get("residual")
    t = 0.02
    bc.integrateResidual(residual, t, fields)

    # No testing of result.
    return


  def test_integrateJacobian(self):
    """
    Test integrateJacobian().

    WARNING: This is not a rigorous test of integrateJacobian() because we
    don't verify the results.
    """

    (mesh, bc, fields) = self._initialize()

    from pylith.topology.Jacobian import Jacobian
    jacobian = Jacobian(fields.solution())
    jacobian.zero()
    t = 0.24
    bc.integrateJacobian(jacobian, t, fields)
    self.assertEqual(False, bc.needNewJacobian())

    # No testing of result.
    return


  def test_poststep(self):
    """
    Test poststep().

    WARNING: This is not a rigorous test of poststep() because we
    neither set the input fields or verify the results.
    """
    (mesh, bc, fields) = self._initialize()

    t = 0.50
    dt = 0.1
    totalTime = 5
    bc.poststep(t, dt, fields)

    # No testing of result.
    return
  

  def test_finalize(self):
    """
    Test finalize().

    WARNING: This is not a rigorous test of finalize() because we
    neither set the input fields or verify the results.
    """
    (mesh, bc, fields) = self._initialize()
    bc.finalize()

    # No testing of result.
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.bc.AbsorbingDampers import boundary_condition
    bc = boundary_condition()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _initialize(self):
    """
    Initialize AbsorbingDampers boundary condition.
    """
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    db = SimpleDB()
    db.inventory.label = "TestAbsorbingDampers tri3"
    db.inventory.iohandler.inventory.filename = \
        "data/elasticplanestrain.spatialdb"
    db.inventory.iohandler._configure()
    db._configure()

    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATSimplex()
    cell.inventory.shape = "line"
    cell.inventory.degree = 1
    cell.inventory.order = 1
    cell._configure()
    from pylith.feassemble.Quadrature import SubMeshQuadrature
    quadrature = SubMeshQuadrature()
    quadrature.inventory.cell = cell
    quadrature._configure()

    from pylith.bc.AbsorbingDampers import AbsorbingDampers
    bc = AbsorbingDampers()
    bc.inventory.quadrature = quadrature
    bc.inventory.db = db
    bc.inventory.id = 0
    bc.inventory.label = "bc"
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
    bc.timeStep(0.01)

    # Setup fields
    from pylith.topology.SolutionFields import SolutionFields
    fields = SolutionFields(mesh)
    fields.add("residual", "residual")
    fields.add("dispIncr(t->t+dt)", "displacement")
    fields.add("disp(t)", "displacement")
    fields.add("disp(t-dt)", "displacement")
    fields.add("velocity(t)", "velocity")
    fields.solutionName("dispIncr(t->t+dt)")

    residual = fields.get("residual")
    residual.newSection(residual.VERTICES_FIELD, cs.spaceDim())
    residual.allocate()
    residual.zero()

    fields.copyLayout("residual")
    
    return (mesh, bc, fields)


# End of file 
