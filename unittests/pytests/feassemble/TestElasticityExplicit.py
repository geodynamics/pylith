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

## @file unittests/pytests/feassemble/TestElasticityExplicit.py

## @brief Unit testing of Python ElasticityExplicit object.

import unittest
from pylith.feassemble.ElasticityExplicit import ElasticityExplicit

from spatialdata.geocoords.CSCart import CSCart

# ----------------------------------------------------------------------
class TestElasticityExplicit(unittest.TestCase):
  """
  Unit testing of Python ElasticityExplicit object.
  """

  def test_implementsIntegrator(self):
    """
    Test to make sure ElasticityExplicit satisfies integrator requirements.
    """
    integrator = ElasticityExplicit()
    from pylith.feassemble.Integrator import implementsIntegrator
    self.failUnless(implementsIntegrator(integrator))
    return
    

  def test_preinitialize(self):
    """
    Test preiniitlaize().
    """
    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    normalizer.initialize()

    # Setup mesh
    cs = CSCart()
    cs.spaceDim = 2
    from pylith.meshio.MeshIOAscii import MeshIOAscii
    importer = MeshIOAscii()
    importer.filename = "data/tri3.mesh"
    importer.coordsys = cs
    mesh = importer.read(normalizer, debug=False, interpolate=False)

    # Setup material
    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATSimplex()
    cell.shape = "triangle"
    cell.degree = 1
    cell.order = 1
    from pylith.feassemble.quadrature.Quadrature2D import Quadrature2D
    quadrature = Quadrature2D()
    quadrature._configure()
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
    initialStateDB = None

    from pylith.materials.ElasticPlaneStrain import ElasticPlaneStrain
    material = ElasticPlaneStrain()
    material.id = 0
    material.label = "elastic plane strain"
    material.db = db
    material.quadrature = quadrature
    material.initialStateDB = initialStateDB
    from pylith.meshio.OutputMatElastic import OutputMatElastic
    material.output = OutputMatElastic()
    material.output._configure()
    material.output.writer._configure()

    integrator = ElasticityExplicit()
    integrator.preinitialize(mesh, material)
    self.assertEqual(mesh, integrator.mesh)
    self.assertEqual(minJacobian, integrator.quadrature.minJacobian)
    return
    

  def test_timeStep(self):
    """
    Test timeStep().
    """
    dt = 2.3
    (mesh, integrator, fields) = self._initialize()
    integrator.timeStep(dt)
    return

  
  def test_stableTimeStep(self):
    """
    Test stableTimeStep().
    """
    (mesh, integrator, fields) = self._initialize()

    self.assertEqual(1.0e+30, integrator.stableTimeStep())
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
    try:
      integrator.useSolnIncr(True)
      self.failIf(True)
    except:
      self.failIf(False)
    return


  def test_integrateResidual(self):
    """
    Test integrateResidual().

    WARNING: This is not a rigorous test of integrateResidual() because we
    neither set the input fields or verify the results.
    """
    (mesh, integrator, fields) = self._initialize()

    residual = fields.getReal("residual")

    t = 0.45
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
    t = 0.145
    integrator.integrateJacobian(jacobian, t, fields)
    self.assertEqual(False, integrator.needNewJacobian())

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return

  
  def test_poststep(self):
    """
    Test poststep().

    WARNING: This is not a rigorous test of poststep() because we
    neither set the input fields or verify the results.
    """
    (mesh, integrator, fields) = self._initialize()

    t = 3.45

    residual = fields.getReal("residual")
    integrator.integrateResidual(residual, t, fields)

    dt = 0.02
    totalTime = 5.0
    integrator.poststep(t, dt, totalTime, fields)

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
    dt = 2.3
    
    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    normalizer.initialize()

    # Setup mesh
    cs = CSCart()
    cs.spaceDim = 2
    from pylith.meshio.MeshIOAscii import MeshIOAscii
    importer = MeshIOAscii()
    importer.filename = "data/tri3.mesh"
    importer.coordsys = cs
    mesh = importer.read(normalizer, debug=False, interpolate=False)

    # Setup material
    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATSimplex()
    cell.shape = "triangle"
    cell.degree = 1
    cell.order = 1
    from pylith.feassemble.quadrature.Quadrature2D import Quadrature2D
    quadrature = Quadrature2D()
    quadrature._configure()
    quadrature.cell = cell
    
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    iohandler = SimpleIOAscii()
    iohandler.filename = "data/elasticplanestrain.spatialdb"
    db = SimpleDB()
    db.label = "elastic plane strain"
    db.iohandler = iohandler
    initialStateDB = None

    from pylith.materials.ElasticPlaneStrain import ElasticPlaneStrain
    material = ElasticPlaneStrain()
    material.id = 0
    material.label = "elastic plane strain"
    material.db = db
    material.initialStateDB = initialStateDB
    material.quadrature = quadrature
    from pylith.meshio.OutputMatElastic import OutputMatElastic
    material.output = OutputMatElastic()
    material.output._configure()
    material.output.writer._configure()

    # Setup integrator
    integrator = ElasticityExplicit()
    integrator.preinitialize(mesh, material)
    from pyre.units.time import s
    integrator.initialize(totalTime=0.0*s, numTimeSteps=1, normalizer=normalizer)
    integrator.timeStep(dt)

    # Setup fields
    from pylith.topology.FieldsManager import FieldsManager
    fields = FieldsManager(mesh)
    fields.addReal("residual")
    fields.addReal("solution")
    fields.addReal("dispT")
    fields.addReal("dispTmdt")
    fields.createHistory(["solution", "dispT", "dispTmdt"])
    fields.solutionField("solution")
    fields.setFiberDimension("residual", cs.spaceDim)
    fields.allocate("residual")
    fields.copyLayout("residual")

    import pylith.topology.topology as bindings
    bindings.zeroRealSection(fields.getReal("residual"))
    
    return (mesh, integrator, fields)


# End of file 
