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

## @file unittests/pytests/bc/TestNeumann.py

## @brief Unit testing of Neumann object.

import unittest

from pylith.bc.Neumann import Neumann

# ----------------------------------------------------------------------
class TestNeumann(unittest.TestCase):
  """
  Unit testing of Neumann object.
  """

  def test_implementsIntegrator(self):
    """
    Test to make sure Neumann satisfies constraint requirements.
    """
    bc = Neumann()
    from pylith.feassemble.Integrator import implementsIntegrator
    self.failUnless(implementsIntegrator(bc))
    return
    

  def test_constructor(self):
    """
    Test constructor.
    """
    from pylith.bc.Neumann import Neumann
    bc = Neumann()
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

    dt = 0.25
    bc.timeStep(dt)
    return

  
  def test_stableTimeStep(self):
    """
    Test stableTimeStep().
    """
    (mesh, bc, fields) = self._initialize()

    self.assertEqual(1.0e+30, bc.stableTimeStep())
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
    t = 0.02
    bc.integrateResidual(residual, t, fields)

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return


  def test_integrateJacobian(self):
    """
    Test integrateJacobian().

    This does nothing for Neumann BC.
    """

    (mesh, bc, fields) = self._initialize()

    jacobian = mesh.createMatrix(fields.getReal("residual"))
    import pylith.utils.petsc as petsc
    petsc.mat_setzero(jacobian)
    t = 0.24
    bc.integrateJacobian(jacobian, t, fields)

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
    bc.poststep(t, dt, totalTime, fields)

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
    Initialize Neumann boundary condition.
    """
    from pylith.bc.Neumann import Neumann
    bc = Neumann()
    bc._configure()
    bc.output._configure()
    bc.output.writer._configure()
    bc.label = "bc"

    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATSimplex()
    cell.shape = "line"
    cell.degree = 1
    cell.order = 1
    from pylith.feassemble.quadrature.Quadrature1Din2D import Quadrature1Din2D
    quadrature = Quadrature1Din2D()
    quadrature._configure()
    quadrature.cell = cell
    bc.quadrature = quadrature

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    db = SimpleDB()
    db._configure()
    db.label = "TestNeumann tri3"
    db.iohandler.filename = "data/tri3-tractions.spatialdb"
    db.initialize()
    bc.db = db

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
    bc.timeStep(0.01)

    # Setup fields
    from pylith.topology.FieldsManager import FieldsManager
    fields = FieldsManager(mesh)
    fields.addReal("residual")
    fields.addReal("dispTBctpdt")
    fields.addReal("dispIncr")
    fields.setFiberDimension("residual", cs.spaceDim)
    fields.allocate("residual")
    fields.copyLayout("residual")

    import pylith.topology.topology as bindings
    bindings.zeroRealSection(fields.getReal("dispIncr"))
    
    return (mesh, bc, fields)


# End of file 
