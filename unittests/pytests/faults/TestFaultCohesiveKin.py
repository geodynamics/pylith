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

  def test_implementsIntegrator(self):
    """
    Test to make sure FaultCohesiveKin satisfies integrator requirements.
    """
    fault = FaultCohesiveKin()
    from pylith.feassemble.Integrator import implementsIntegrator
    self.assertTrue(implementsIntegrator(fault))
    return
    

  def test_constructor(self):
    """
    Test constructor.
    """
    fault = FaultCohesiveKin()
    self.failIfEqual(None, fault.cppHandle)
    return


  def test_adjustTopology(self):
    """
    Test adjustTopology().

    WARNING: This is not a rigorous test of updateState() because we
    neither set the input fields or verify the results.
    """
    cs = CSCart()
    cs.spaceDim = 2
    
    from pylith.meshio.MeshIOAscii import MeshIOAscii
    importer = MeshIOAscii()
    importer.filename = "data/tri3.mesh"
    importer.coordsys = cs
    mesh = importer.read(debug=False, interpolate=False)

    fault = FaultCohesiveKin()
    fault.id = 10
    fault.label = "fault"

    fault.adjustTopology(mesh)

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return


  def test_initialize(self):
    """
    Test initialize().

    WARNING: This is not a rigorous test of updateState() because we
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
    from pyre.units.time import second
    dt = 2.4*second
    fault = FaultCohesiveKin()
    fault.timeStep = dt
    self.assertEqual(dt, fault.timeStep)
    return

  
  def test_stableTimeStep(self):
    """
    Test stableTimeStep().
    """
    from pyre.units.time import second
    dt = 2.3*second
    fault = FaultCohesiveKin()
    fault.timeStep(dt)
    self.assertEqual(dt, fault.stableTimeStep())
    return

  
  def test_needNewJacobian(self):
    """
    Test needNewJacobian().
    """
    fault = FaultCohesiveKin()
    self.assertEqual(False, fault.needNewJacobian())
    return

  
  def test_integrateResidual(self):
    """
    Test integrateResidual().

    WARNING: This is not a rigorous test of integrateResidual() because we
    neither set the input fields or verify the results.
    """
    (mesh, fault, fields) = self._initialize()

    residual = fields.getReal("residual")
    t = 1.0*second
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

    jacobian = mesh.createMatrix(fields.getReal("residual"))
    import pylith.utils.petsc as petsc
    petsc.mat_setzero(jacobian)
    t = 1.0*second
    fault.integrateJacobian(jacobian, t, fields)

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return

  
  def test_updateState(self):
    """
    Test updateState().

    WARNING: This is not a rigorous test of updateState() because we
    neither set the input fields or verify the results.
    """
    (mesh, fault, fields) = self._initialize()

    disp = fields.getReal("dispT")
    t = 1.0*second
    fault.updateState(t, disp)

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
    from pyre.units.time import second
    dt = 1.0*second
    
    # Setup mesh
    cs = CSCart()
    cs.spaceDim = 2
    from pylith.meshio.MeshIOAscii import MeshIOAscii
    importer = MeshIOAscii()
    importer.filename = "data/tri3.mesh"
    importer.coordsys = cs
    mesh = importer.read(debug=False, interpolate=False)

    # Setup quadrature
    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATSimplex()
    cell.shape = "line"
    cell.degree = 1
    cell.order = 1
    from pylith.feassemble.quadrature.Quadrature1Din2D import Quadrature1Din2D
    quadrature = Quadrature1Din2D()
    quadrature.cell = cell

    # Setup earthquake source
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    ioFinalSlip = SimpleIOAscii()
    ioFinalSlip.filename = "data/tri3_finalslip.spatialdb"
    dbFinalSlip = SimpleDB()
    dbFinalSlip.iohandler = ioFinalSlip
    dbFinalSlip.label = "final slip"
    
    ioSlipTime = SimpleIOAscii()
    ioSlipTime.filename = "data/tri3_sliptime.spatialdb"
    dbSlipTime = SimpleDB()
    dbSlipTime.iohandler = ioSlipTime
    dbSlipTime.label = "slip time"
    
    ioPeakRate = SimpleIOAscii()
    ioPeakRate.filename = "data/tri3_peakrate.spatialdb"
    dbPeakRate = SimpleDB()
    dbPeakRate.iohandler = ioPeakRate
    dbPeakRate.label = "peak rate"
    
    from pylith.faults.BruneSlipFn import BruneSlipFn
    slipfn = BruneSlipFn()
    slipfn.slip = dbFinalSlip
    slipfn.slipTime = dbSlipTime
    slipfn.slipRate = dbPeakRate

    from pylith.faults.EqKinSrc import EqKinSrc
    eqsrc = EqKinSrc()
    eqsrc.slipfn = slipfn

    # Setup fault
    fault = FaultCohesiveKin()
    fault.id = 10
    fault.label = "fault"
    fault.upDir = [0, 0, 1]
    fault.quadrature = quadrature
    fault.eqsrc = eqsrc
    fault.timeStep(dt)
    fault.adjustTopology(mesh)
    fault.initialize(mesh)

    # Setup fields
    from pylith.topology.FieldsManager import FieldsManager
    fields = FieldsManager(mesh)
    fields.addReal("residual")
    fields.addReal("dispT")
    fields.solutionField("dispT");
    fields.setFiberDimension("residual", cs.spaceDim)
    fields.allocate("residual")
    fields.copyLayout("residual")

    import pylith.topology.topology as bindings
    bindings.zeroRealSection(fields.getReal("residual"))
    
    return (mesh, fault, fields)


# End of file 
