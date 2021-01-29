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

## @file tests/pytests/faults/TestFaultCohesiveDyn.py

## @brief Unit testing of FaultCohesiveDyn object.

import unittest

from pylith.faults.FaultCohesiveDyn import FaultCohesiveDyn

from spatialdata.geocoords.CSCart import CSCart
from pythia.pyre.units.time import second

# ----------------------------------------------------------------------
class TestFaultCohesiveDyn(unittest.TestCase):
  """
  Unit testing of Fault object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    fault = FaultCohesiveDyn()
    return


  def test_configure(self):
    """
    Test _configure().
    """
    fault = FaultCohesiveDyn()
    fault.inventory.faultLabel = "fault group"
    fault._configure()
    return


  def test_implementsIntegrator(self):
    """
    Test to make sure FaultCohesiveDyn satisfies integrator requirements.
    """
    fault = FaultCohesiveDyn()
    from pylith.feassemble.Integrator import implementsIntegrator
    self.failUnless(implementsIntegrator(fault))
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

    fault = FaultCohesiveDyn()
    fault.inventory.matId = 10
    fault.inventory.faultLabel = "fault"
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
    from pylith.faults.FaultCohesiveDyn import fault
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

    # Setup rupture info
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    ioTractions = SimpleIOAscii()
    ioTractions.inventory.filename = "data/tri3_initialtractions.spatialdb"
    ioTractions._configure()
    dbTractions = SimpleDB()
    dbTractions.inventory.iohandler = ioTractions
    dbTractions.inventory.label = "initial tractions"
    dbTractions._configure()
    from pylith.faults.TractPerturbation import TractPerturbation
    tract = TractPerturbation()
    tract.inventory.dbInitial = dbTractions
    tract._configure()

    ioFriction = SimpleIOAscii()
    ioFriction.inventory.filename = "data/tri3_staticfriction.spatialdb"
    ioFriction._configure()
    dbFriction = SimpleDB()
    dbFriction.inventory.iohandler = ioFriction
    dbFriction.inventory.label = "friction"
    dbFriction._configure()
    
    from pylith.friction.StaticFriction import StaticFriction
    friction = StaticFriction()
    friction.inventory.label = "Static friction"
    friction.inventory.dbProperties = dbFriction
    friction._configure()

    # Setup fault
    fault = FaultCohesiveDyn()
    fault.inventory.output.inventory.writer._configure()
    fault.inventory.output._configure()
    fault.inventory.matId = 10
    fault.inventory.faultLabel = "fault"
    fault.inventory.upDir = [0, 0, 1]
    fault.inventory.faultQuadrature = quadrature
    fault.inventory.tract = tract
    fault.inventory.friction = friction
    fault._configure()

    nvertices = fault.numVerticesNoMesh(mesh)
    firstFaultVertex = 0
    firstLagrangeVertex = nvertices
    firstFaultCell      = 2*nvertices
    fault.adjustTopology(mesh, firstFaultVertex, firstLagrangeVertex,
                         firstFaultCell)
    from pylith.topology.topology import MeshOps_nondimensionalize
    MeshOps_nondimensionalize(mesh, normalizer)

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
    fields.add("velocity(t)", "velocity")
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
