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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file unittests/pytests/feassemble/TestElasticityExplicitLgDeform.py

## @brief Unit testing of Python ElasticityExplicitLgDeform object.

import unittest
from pylith.feassemble.ElasticityExplicitLgDeform import ElasticityExplicitLgDeform

from spatialdata.geocoords.CSCart import CSCart

# ----------------------------------------------------------------------
class TestElasticityExplicitLgDeform(unittest.TestCase):
  """
  Unit testing of Python ElasticityExplicitLgDeform object.
  """

  def test_implementsIntegrator(self):
    """
    Test to make sure ElasticityExplicitLgDeform satisfies integrator requirements.
    """
    integrator = ElasticityExplicitLgDeform()
    from pylith.feassemble.Integrator import implementsIntegrator
    self.failUnless(implementsIntegrator(integrator))
    return
    

  def test_preinitialize(self):
    """
    Test preiniitlaize().

    WARNING: This is not a rigorous test of preinitialize() because we
    neither set the input fields or verify the results.
    """
    (mesh, integrator) = self._preinitialize()

    # No test of result.
    return


  def test_verifyConfiguration(self):
    """
    Test verifyConfiguration().

    WARNING: This is not a rigorous test of verifyConfiguration()
    because we neither set the input fields or verify the results.
    """
    (mesh, integrator) = self._preinitialize()
    integrator.verifyConfiguration()

    # No test of result.
    return


  def test_initialize(self):
    """
    Test initialize().

    WARNING: This is not a rigorous test of initialize() because we
    neither set the input fields or verify the results.
    """
    (mesh, integrator) = self._preinitialize()
    fields = self._initialize(mesh, integrator)

    # No test of result.
    return


  def test_timeStep(self):
    """
    Test timeStep().

    WARNING: This is not a rigorous test of timeStep() because we
    neither set the input fields or verify the results.
    """
    dt = 2.3
    (mesh, integrator) = self._preinitialize()
    fields = self._initialize(mesh, integrator)
    integrator.timeStep(dt)

    # No test of result.
    return


  def test_normViscosity(self):
    """
    Test normViscosity().

    WARNING: This is not a rigorous test of timeStep() because we
    neither set the input fields or verify the results.
    """
    viscosity = 1.234
    (mesh, integrator) = self._preinitialize()
    fields = self._initialize(mesh, integrator)
    integrator.normViscosity(viscosity)

    # No test of result.
    return


  def test_stableTimeStep(self):
    """
    Test stableTimeStep().
    """
    (mesh, integrator) = self._preinitialize()
    fields = self._initialize(mesh, integrator)

    from pylith.utils.utils import maxscalar
    self.assertAlmostEqual(1.0, integrator.stableTimeStep(mesh)/maxscalar(), 7)
    return


  def test_needNewJacobian(self):
    """
    Test needNewJacobian().
    """
    (mesh, integrator) = self._preinitialize()
    fields = self._initialize(mesh, integrator)

    self.assertEqual(True, integrator.needNewJacobian())
    return


  def test_useSolnIncr(self):
    """
    Test useSolnIncr().

    WARNING: This is not a rigorous test of useSolnIncr() because we
    neither set the input fields or verify the results.
    """
    (mesh, integrator) = self._preinitialize()
    fields = self._initialize(mesh, integrator)

    try:
      integrator.useSolnIncr(True)
      self.failIf(True)
    except:
      self.failIf(False)

    # No test of result.
    return


  def test_integrateResidual(self):
    """
    Test integrateResidual().

    WARNING: This is not a rigorous test of integrateResidual()
    because we neither set the input fields or verify the results.
    """
    (mesh, integrator) = self._preinitialize()
    fields = self._initialize(mesh, integrator)

    residual = fields.get("residual")
    t = 3.4
    integrator.integrateResidual(residual, t, fields)

    # No test of result.
    return


  def test_integrateJacobian(self):
    """
    Test integrateJacobian().

    WARNING: This is not a rigorous test of integrateJacobian()
    because we neither set the input fields or verify the results.
    """
    (mesh, integrator) = self._preinitialize()
    fields = self._initialize(mesh, integrator)

    from pylith.topology.Field import Field
    jacobian = Field(mesh)
    spaceDim = mesh.coordsys().spaceDim()
    jacobian.subfieldAdd("displacement", spaceDim, jacobian.VECTOR);
    jacobian.subfieldAdd("lagrange_multiplier", spaceDim, jacobian.VECTOR);

    jacobian.subfieldsSetup();
    jacobian.setupSolnChart();
    jacobian.setupSolnDof(spaceDim);
    jacobian.allocate();
    jacobian.zeroAll();

    t = 7.3
    self.assertEqual(True, integrator.needNewJacobian())
    integrator.integrateJacobian(jacobian, t, fields)
    self.assertEqual(False, integrator.needNewJacobian())
    
    # No test of result.
    return


  def test_poststep(self):
    """
    Test poststep().

    WARNING: This is not a rigorous test of poststep() because we
    neither set the input fields or verify the results.
    """
    (mesh, integrator) = self._preinitialize()
    fields = self._initialize(mesh, integrator)

    t = 7.3
    dt = 0.1
    totalTime = 23.0
    integrator.poststep(t, dt, fields)

    # No test of result
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.feassemble.ElasticityExplicitLgDeform import integrator
    i = integrator()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _preinitialize(self):
    """
    Setup mesh and integrator and preinitialize integrator.
    """
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

    # Setup material
    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATSimplex()
    cell.inventory.dimension = 2
    cell.inventory.degree = 1
    cell.inventory.order = 1
    cell._configure()
    from pylith.feassemble.Quadrature import Quadrature
    quadrature = Quadrature()
    quadrature.inventory.cell = cell
    quadrature._configure()
    
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    iohandler = SimpleIOAscii()
    iohandler.inventory.filename = "data/elasticplanestrain.spatialdb"
    iohandler._configure()
    db = SimpleDB()
    db.inventory.label = "elastic plane strain"
    db.inventory.iohandler = iohandler
    db._configure()

    from pylith.materials.ElasticPlaneStrain import ElasticPlaneStrain
    material = ElasticPlaneStrain()
    material.inventory.label = "elastic plane strain"
    material.inventory.id = 0
    material.inventory.dbProperties = db
    material.inventory.quadrature = quadrature
    material._configure()
    
    from pylith.meshio.OutputMatElastic import OutputMatElastic
    material.output = OutputMatElastic()
    material.output._configure()
    material.output.writer._configure()

    # Setup integrator
    integrator = ElasticityExplicitLgDeform()
    integrator.preinitialize(mesh, material)
    return (mesh, integrator)


  def _initialize(self, mesh, integrator):
    """
    Initialize integrator.
    """
    dt = 2.3
    
    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    normalizer._configure()

    from pyre.units.time import s
    integrator.initialize(totalTime=0.0*s, numTimeSteps=1,
                          normalizer=normalizer)
    integrator.timeStep(dt)

    # Setup fields
    from pylith.topology.SolutionFields import SolutionFields
    fields = SolutionFields(mesh)
    fields.add("residual", "residual")
    fields.add("disp(t+dt)", "displacement")
    fields.add("disp(t)", "displacement")
    fields.add("disp(t-dt)", "displacement")
    fields.add("velocity(t)", "velocity")
    fields.add("acceleration(t)", "acceleration")
    fields.solutionName("disp(t+dt)")

    residual = fields.get("residual")
    spaceDim = mesh.coordsys().spaceDim()
    lengthScale = normalizer.lengthScale()
    residual.subfieldAdd("displacement", spaceDim, residual.VECTOR, lengthScale.value);
    residual.subfieldAdd("lagrange_multiplier", spaceDim, residual.VECTOR);
    
    residual.subfieldsSetup();
    residual.setupSolnChart();
    residual.setupSolnDof(spaceDim);
    residual.allocate();
    residual.zeroAll();
    fields.copyLayout("residual")

    return fields


# End of file 
