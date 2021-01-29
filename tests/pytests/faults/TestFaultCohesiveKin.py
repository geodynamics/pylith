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
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file tests/pytests/faults/TestFaultCohesiveKin.py

## @brief Unit testing of FaultCohesiveKin object.

import unittest

from pylith.faults.FaultCohesiveKin import FaultCohesiveKin

from spatialdata.geocoords.CSCart import CSCart
from pythia.pyre.units.time import second

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
    fault.inventory.faultLabel = "fault group"
    fault.inventory.faultEdge = "fault edge"
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
    fault.inventory.faultLabel = "fault group"
    fault._configure()

    fault.useFaultMesh(True);

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
    mesh = importer.read(debug=False, interpolate=False)

    fault = FaultCohesiveKin()
    fault.inventory.matId = 10
    fault.inventory.faultLabel = "fault"
    fault.inventory.faultEdge = "fault_edge"
    fault._configure()

    nvertices = fault.numVerticesNoMesh(mesh)
    firstFaultVertex = 0
    firstLagrangeVertex = nvertices
    firstFaultCell      = 2*nvertices
    fault.adjustTopology(mesh, firstFaultVertex, firstLagrangeVertex,
                         firstFaultCell)

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

    from pylith.utils.utils import maxscalar
    self.assertAlmostEqual(1.0, fault.stableTimeStep(mesh)/maxscalar(), 7)
    return

  
  def test_needNewJacobian(self):
    """
    Test needNewJacobian().
    """
    (mesh, fault, fields) = self._initialize()
    self.assertEqual(True, fault.needNewJacobian())
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
    jacobian = Jacobian(fields.solution())
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
    fault.poststep(t, dt, fields)

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


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.faults.FaultCohesiveKin import fault
    f = fault()
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
    mesh = importer.read(debug=False, interpolate=False)

    # Setup quadrature
    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATSimplex()
    cell.inventory.dimension = 1
    cell.inventory.degree = 1
    cell.inventory.order = 1
    cell._configure()
    from pylith.feassemble.Quadrature import Quadrature
    quadrature = Quadrature()
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

    # Setup fault
    fault = FaultCohesiveKin()
    fault.inventory.output.inventory.writer._configure()
    fault.inventory.output._configure()
    fault.inventory.matId = 10
    fault.inventory.faultLabel = "fault"
    fault.inventory.upDir = [0, 0, 1]
    fault.inventory.faultQuadrature = quadrature
    fault._configure()
    eqsrc = fault.eqsrcs.components()[0]
    eqsrc.inventory.originTime = 1.23*second
    eqsrc.inventory.slipfn = slipfn
    eqsrc._configure()

    nvertices = fault.numVerticesNoMesh(mesh)
    firstFaultVertex = 0
    firstLagrangeVertex = nvertices
    firstFaultCell      = 2*nvertices
    fault.adjustTopology(mesh, firstFaultVertex, firstLagrangeVertex,
                         firstFaultCell)
    fault.preinitialize(mesh)
    fault.timeStep(dt)
    fault.verifyConfiguration()
    from pythia.pyre.units.time import s
    fault.initialize(totalTime=0.0*s, numTimeSteps=1, normalizer=normalizer)

    # Setup fields
    from pylith.topology.SolutionFields import SolutionFields
    fields = SolutionFields(mesh)
    fields.add("residual", "residual")
    fields.add("dispIncr(t->t+dt)", "displacement_increment")
    fields.add("disp(t)", "displacement")
    fields.solutionName("dispIncr(t->t+dt)")

    residual = fields.get("residual")
    residual.subfieldAdd("displacement", cs.getSpaceDim(), residual.VECTOR)
    residual.subfieldAdd("lagrange_multiplier", cs.getSpaceDim(), residual.VECTOR)
    residual.subfieldsSetup()
    residual.setupSolnChart()
    residual.setupSolnDof(cs.getSpaceDim())
    fault.setupSolnDof(residual)
    residual.allocate()
    residual.zero()

    fields.copyLayout("residual")
    
    return (mesh, fault, fields)


# End of file 
