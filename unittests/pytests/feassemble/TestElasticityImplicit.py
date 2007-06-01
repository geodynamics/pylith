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

# ----------------------------------------------------------------------
class TestElasticityImplicit(unittest.TestCase):
  """
  Unit testing of Python ElasticityImplicit object.
  """

  def test_initQuadrature(self):
    """
    Test initQuadrature().
    """
    from pylith.feassemble.quadrature.Quadrature2D import Quadrature2D
    q = Quadrature2D()
    minJacobian = 4.0e-02;
    q.minJacobian = minJacobian
    from pylith.feassemble.FIATSimplex import FIATSimplex
    q.cell = FIATSimplex()
    q.cell.shape = "triangle"
    q.cell.order = 1
    q.cell.degree = 1

    integrator = ElasticityImplicit()
    integrator.initQuadrature(q)
    self.assertEqual(minJacobian, integrator.quadrature.minJacobian)
    return
    

  def test_needNewJacobian(self):
    """
    Test needNewJacobian().
    """
    integrator = ElasticityImplicit()
    self.assertEqual(False, integrator.needNewJacobian())
    return


  def test_setMesh(self):
    """
    Test setMesh().
    """
    cs = CSCart()
    cs.spaceDim = 2
    
    from pylith.meshio.MeshIOAscii import MeshIOAscii
    importer = MeshIOAscii()
    importer.filename = "data/tri3.mesh"
    importer.coordsys = cs
    mesh = importer.read(debug=False, interpolate=False)

    integrator = ElasticityImplicit()
    integrator.setMesh(mesh)
    self.assertEqual(mesh, integrator.mesh)
    return

  
  def test_timeStep(self):
    """
    Test timeStep().
    """
    from pyre.units.time import second
    dt = 2.4*second
    integrator = ElasticityImplicit()
    integrator.timeStep = dt
    self.assertEqual(dt, integrator.timeStep)
    return

  
  def test_stableTimeStep(self):
    """
    Test stableTimeStep().
    """
    from pyre.units.time import second
    dt = 2.3*second
    integrator = ElasticityImplicit()
    integrator.timeStep(dt)
    self.assertEqual(dt, integrator.stableTimeStep())
    return

  
  def test_needNewJacobian(self):
    """
    Test needNewJacobian().
    """
    integrator = ElasticityImplicit()
    self.assertEqual(False, integrator.needNewJacobian())
    return

  
  def test_integrateResidual(self):
    """
    Test integrateResidual().

    WARNING: This is not a rigorous test of integrateResidual() because we
    neither set the input fields or verify the results.
    """
    (mesh, integrator, fields) = self._initialize()

    residual = fields.getReal("residual")
    integrator.integrateResidual(residual, fields)

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
    integrator.integrateJacobian(jacobian, fields)

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

    dispT = fields.getReal("dispT")
    integrator.updateState(dispT)

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _initialize(self):
    """
    Initialize integrator.
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
    integrator.setMesh(mesh)
    integrator.initQuadrature(material.quadrature)
    integrator.initMaterial(mesh, material)
    integrator.timeStep(dt)

    # Setup fields
    from pylith.topology.FieldsManager import FieldsManager
    fields = FieldsManager(mesh)
    fields.addReal("residual")
    fields.addReal("dispBCTpdt")
    fields.addReal("dispT")
    fields.createHistory(["dispBCTpdt", "dispT"])
    fields.setFiberDimension("residual", cs.spaceDim)
    fields.allocate("residual")
    fields.copyLayout("residual")

    import pylith.topology.topology as bindings
    bindings.zeroRealSection(fields.getReal("residual"))
    
    return (mesh, integrator, fields)


# End of file 
