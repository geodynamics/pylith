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

## @file unittests/pytests/faults/TestFaultCohesiveKin.py

## @brief Unit testing of FaultCohesiveKin object.

import unittest

from pylith.faults.FaultCohesiveKin import FaultCohesiveKin

from spatialdata.geocoords.CSCart import CSCart
from pyre.units.time import second

# ----------------------------------------------------------------------
class TestFaultCohesiveKin(unittest.TestCase):
  """
  Unit testing of Fault object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    fault = FaultCohesiveKin()
    return


  def test_configure(self):
    """
    Test _configure().
    """
    fault = FaultCohesiveKin()
    fault._configure()
    return


  def test_implementsIntegrator(self):
    """
    Test to make sure FaultCohesiveKin satisfies integrator requirements.
    """
    fault = FaultCohesiveKin()
    from pylith.feassemble.Integrator import implementsIntegrator
    self.failUnless(implementsIntegrator(fault))
    return
    

  def test_useFaultMesh(self):
    """
    Test useFaultMesh().
    """
    fault = FaultCohesiveKin()
    fault._configure()

    fault.useFaultMesh(True);

    # No test of result
    return


  def test_faultMeshFilename(self):
    """
    Test faultMeshFilename().
    """
    fault = FaultCohesiveKin()
    fault._configure()

    filename = "SanAndreas.inp"
    fault.faultMeshFilename(filename)

    # No test of result
    return


  def test_adjustTopology(self):
    """
    Test adjustTopology().

    WARNING: This is not a rigorous test of adjustTopology() because we
    neither set the input fields or verify the results.
    """
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

    fault = FaultCohesiveKin()
    fault.inventory.matId = 10
    fault.inventory.faultLabel = "fault"
    fault._configure()

    fault.adjustTopology(mesh)

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return


  def test_initialize(self):
    """
    Test initialize().

    WARNING: This is not a rigorous test of initialize() because we
    neither set the input fields or verify the results.
    """
    (mesh, fault, fields) = self._initialize()
    
    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return


  def test_timeStep(self):
    """
    Test timeStep().
    """
    dt = 2.4
    (mesh, fault, fields) = self._initialize()
    fault.timeStep(dt)

    # No test of result
    return

  
  def test_stableTimeStep(self):
    """
    Test stableTimeStep().
    """
    (mesh, fault, fields) = self._initialize()

    self.assertEqual(1.0e+30, fault.stableTimeStep(mesh))
    return

  
  def test_needNewJacobian(self):
    """
    Test needNewJacobian().
    """
    (mesh, fault, fields) = self._initialize()
    self.assertEqual(True, fault.needNewJacobian())
    return

  
  def test_useSolnIncr(self):
    """
    Test useSolnIncr().
    """
    (mesh, fault, fields) = self._initialize()
    fault.useSolnIncr(True)

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return

  
  def test_integrateResidual(self):
    """
    Test integrateResidual().

    WARNING: This is not a rigorous test of integrateResidual() because we
    neither set the input fields or verify the results.
    """
    (mesh, fault, fields) = self._initialize()

    residual = fields.get("residual")
    t = 1.0
    fault.integrateResidual(residual, t, fields)

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return

  
  def test_integrateJacobian(self):
    """
    Test integrateJacobian().

    WARNING: This is not a rigorous test of integrateJacobian() because we
    neither set the input fields or verify the results.
    """
    (mesh, fault, fields) = self._initialize()

    from pylith.topology.Jacobian import Jacobian
    jacobian = Jacobian(fields)
    jacobian.zero()
    t = 1.0
    fault.integrateJacobian(jacobian, t, fields)
    self.assertEqual(False, fault.needNewJacobian())

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return

  
  def test_poststep(self):
    """
    Test poststep().

    WARNING: This is not a rigorous test of poststep() because we
    neither set the input fields or verify the results.
    """
    (mesh, fault, fields) = self._initialize()

    t = 0.50
    residual = fields.get("residual")
    fault.integrateResidual(residual, t, fields)

    dt = 0.1
    totalTime = 5
    fault.poststep(t, dt, totalTime, fields)

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return
  

  def test_finalize(self):
    """
    Test finalize().

    WARNING: This is not a rigorous test of finalize() because we
    neither set the input fields or verify the results.
    """
    (mesh, fault, fields) = self._initialize()

    fault.finalize()

    # We should really add something here to check to make sure things
    # actually initialized correctly.
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _initialize(self):
    """
    Initialize fault.
    """
    dt = 2.4
    
    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    normalizer._configure()

    # Setup mesh
    cs = CSCart()
    cs.inventory.spaceDim = 2
    cs._configure()

    from pylith.meshio.MeshIOAscii import MeshIOAscii
    importer = MeshIOAscii()
    importer.inventory.filename = "data/tri3.mesh"
    importer.inventory.coordsys = cs
    importer._configure()
    mesh = importer.read(normalizer, debug=False, interpolate=False)

    # Setup quadrature
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

    # Setup earthquake source
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    ioFinalSlip = SimpleIOAscii()
    ioFinalSlip.inventory.filename = "data/tri3_finalslip.spatialdb"
    ioFinalSlip._configure()
    dbFinalSlip = SimpleDB()
    dbFinalSlip.inventory.iohandler = ioFinalSlip
    dbFinalSlip.inventory.label = "final slip"
    dbFinalSlip._configure()
    
    ioSlipTime = SimpleIOAscii()
    ioSlipTime.inventory.filename = "data/tri3_sliptime.spatialdb"
    ioSlipTime._configure()
    dbSlipTime = SimpleDB()
    dbSlipTime.inventory.iohandler = ioSlipTime
    dbSlipTime.inventory.label = "slip time"
    dbSlipTime._configure()
    
    from pylith.faults.StepSlipFn import StepSlipFn
    slipfn = StepSlipFn()
    slipfn.inventory.dbSlip = dbFinalSlip
    slipfn.inventory.dbSlipTime = dbSlipTime
    slipfn._configure()

    ioMatDB = SimpleIOAscii()
    ioMatDB.inventory.filename = "data/bulkprops_2d.spatialdb"
    ioMatDB._configure()
    dbMat = SimpleDB()
    dbMat.inventory.iohandler = ioMatDB
    dbMat.inventory.label = "bulk properties"
    dbMat._configure()
    
    # Setup fault
    fault = FaultCohesiveKin()
    fault.inventory.output.inventory.writer._configure()
    fault.inventory.output._configure()
    fault.inventory.matId = 10
    fault.inventory.faultLabel = "fault"
    fault.inventory.upDir = [0, 0, 1]
    fault.inventory.normalDir = [1, 0, 0]
    fault.inventory.quadrature = quadrature
    fault.inventory.matDB = dbMat
    fault._configure()
    eqsrc = fault.eqsrcs.components()[0]
    eqsrc.inventory.originTime = 1.23*second
    eqsrc.inventory.slipfn = slipfn
    eqsrc._configure()

    fault.adjustTopology(mesh)
    fault.preinitialize(mesh)
    fault.timeStep(dt)
    fault.verifyConfiguration()
    from pyre.units.time import s
    fault.initialize(totalTime=0.0*s, numTimeSteps=1, normalizer=normalizer)

    # Setup fields
    from pylith.topology.SolutionFields import SolutionFields
    fields = SolutionFields(mesh)
    fields.add("residual", "residual")
    fields.add("solution", "displacement")
    fields.add("disp", "displacement")
    fields.solutionName("solution")
    residual = fields.get("residual")
    residual.newSection(residual.VERTICES_FIELD, cs.spaceDim())
    residual.allocate()
    residual.zero()
    fields.copyLayout("residual")
    
    return (mesh, fault, fields)


# End of file 
