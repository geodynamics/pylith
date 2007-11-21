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

## @file unittests/pytests/feassemble/TestElasticityImplicit.py

## @brief Unit testing of Python ElasticityImplicit object.

import unittest
from pylith.feassemble.ElasticityImplicit import ElasticityImplicit

from spatialdata.geocoords.CSCart import CSCart
from pyre.units.time import year

# ----------------------------------------------------------------------
class TestElasticityImplicit(unittest.TestCase):
  """
  Unit testing of Python ElasticityImplicit object.
  """

  def test_implementsIntegrator(self):
    """
    Test to make sure ElasticityImplicit satisfies integrator requirements.
    """
    integrator = ElasticityImplicit()
    from pylith.feassemble.Integrator import implementsIntegrator
    self.failUnless(implementsIntegrator(integrator))
    return
    

  def test_preinitialize(self):
    """
    Test preiniitlaize().
    """
    # Setup mesh
    cs = CSCart()
    cs.spaceDim = 2
    from pylith.meshio.MeshIOAscii import MeshIOAscii
    importer = MeshIOAscii()
    importer.filename = "data/tri3.mesh"
    importer.coordsys = cs
    mesh = importer.read(debug=False, interpolate=False)

    # Setup material
    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATSimplex()
    cell.shape = "triangle"
    cell.degree = 1
    cell.order = 1
    from pylith.feassemble.quadrature.Quadrature2D import Quadrature2D
    quadrature = Quadrature2D()
    quadrature.cell = cell
    minJacobian = 4.0e-02;
    quadrature.minJacobian = minJacobian
    
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    iohandler = SimpleIOAscii()
    iohandler.filename = "data/elasticplanestrain.spatialdb"
    db = SimpleDB()
    db.label = "elastic plane strain"
    db.iohandler = iohandler

    from pylith.materials.ElasticPlaneStrain import ElasticPlaneStrain
    material = ElasticPlaneStrain()
    material.id = 0
    material.label = "elastic plane strain"
    material.db = db
    material.quadrature = quadrature

    integrator = ElasticityImplicit()
    integrator.preinitialize(mesh, material)
    self.assertEqual(mesh, integrator.mesh)
    self.assertEqual(minJacobian, integrator.quadrature.minJacobian)
    return
    

  def test_timeStep(self):
    """
    Test timeStep().
    """
    from pyre.units.time import second
    dt = 2.3*second
    (mesh, integrator, fields) = self._initialize()
    integrator.timeStep(dt)
    return

  
  def test_stableTimeStep(self):
    """
    Test stableTimeStep().
    """
    from pyre.units.time import second
    dt = 2.3*second
    (mesh, integrator, fields) = self._initialize()
    integrator.timeStep(dt)
    self.assertEqual(dt, integrator.stableTimeStep())
    return

  
  def test_needNewJacobian(self):
    """
    Test needNewJacobian().
    """
    (mesh, integrator, fields) = self._initialize()
    self.assertEqual(True, integrator.needNewJacobian())
    return

  
  def test_useSolnIncr(self):
    """
    Test useSolnIncr().
    """
    (mesh, integrator, fields) = self._initialize()
    integrator.useSolnIncr(True)
    return


  def test_integrateResidual(self):
    """
    Test integrateResidual().

    WARNING: This is not a rigorous test of integrateResidual() because we
    neither set the input fields or verify the results.
    """
    (mesh, integrator, fields) = self._initialize()

    residual = fields.getReal("residual")
    t = 3.4*year
    integrator.integrateResidual(residual, t, fields)

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return

  
  def test_integrateJacobian(self):
    """
    Test integrateJacobian().

    WARNING: This is not a rigorous test of integrateJacobian() because we
    neither set the input fields or verify the results.
    """
    (mesh, integrator, fields) = self._initialize()

    jacobian = mesh.createMatrix(fields.getReal("residual"))
    import pylith.utils.petsc as petsc
    petsc.mat_setzero(jacobian)
    t = 7.3*year
    integrator.integrateJacobian(jacobian, t, fields)
    self.assertEqual(False, integrator.needNewJacobian())

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return

  
  def test_updateState(self):
    """
    Test updateState().

    WARNING: This is not a rigorous test of updateState() because we
    neither set the input fields or verify the results.
    """
    (mesh, integrator, fields) = self._initialize()

    dispTBctpdt = fields.getReal("dispTBctpdt")
    t = 0.27*year
    integrator.updateState(t, dispTBctpdt)

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return
  

  def test_finalize(self):
    """
    Test finalize().

    WARNING: This is not a rigorous test of finalize() because we
    neither set the input fields or verify the results.
    """
    (mesh, integrator, fields) = self._initialize()

    integrator.finalize()

    # We should really add something here to check to make sure things
    # actually initialized correctly.
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _initialize(self):
    """
    Initialize integrator.
    """
    from pyre.units.time import second
    dt = 2.3*second
    
    # Setup mesh
    cs = CSCart()
    cs.spaceDim = 2
    from pylith.meshio.MeshIOAscii import MeshIOAscii
    importer = MeshIOAscii()
    importer.filename = "data/tri3.mesh"
    importer.coordsys = cs
    mesh = importer.read(debug=False, interpolate=False)

    # Setup material
    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATSimplex()
    cell.shape = "triangle"
    cell.degree = 1
    cell.order = 1
    from pylith.feassemble.quadrature.Quadrature2D import Quadrature2D
    quadrature = Quadrature2D()
    quadrature.cell = cell
    
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    iohandler = SimpleIOAscii()
    iohandler.filename = "data/elasticplanestrain.spatialdb"
    db = SimpleDB()
    db.label = "elastic plane strain"
    db.iohandler = iohandler

    from pylith.materials.ElasticPlaneStrain import ElasticPlaneStrain
    material = ElasticPlaneStrain()
    material.id = 0
    material.label = "elastic plane strain"
    material.db = db
    material.quadrature = quadrature

    # Setup integrator
    integrator = ElasticityImplicit()
    integrator.preinitialize(mesh, material)
    integrator.initialize()
    integrator.timeStep(dt)

    # Setup fields
    from pylith.topology.FieldsManager import FieldsManager
    fields = FieldsManager(mesh)
    fields.addReal("residual")
    fields.addReal("dispTBctpdt")
    fields.setFiberDimension("residual", cs.spaceDim)
    fields.allocate("residual")
    fields.copyLayout("residual")

    import pylith.topology.topology as bindings
    bindings.zeroRealSection(fields.getReal("residual"))
    
    return (mesh, integrator, fields)


# End of file 
