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
    self._preinitialize()
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
    mesh = importer.read(normalizer, debug=False, interpolate=False)

    # Setup material
    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATSimplex()
    cell.inventory.shape = "triangle"
    cell.inventory.degree = 1
    cell.inventory.order = 1
    cell._configure()
    from pylith.feassemble.Quadrature import MeshQuadrature
    quadrature = MeshQuadrature()
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
    material.inventory.db = db
    material.inventory.quadrature = quadrature
    material._configure()
    
    from pylith.meshio.OutputMatElastic import OutputMatElastic
    material.output = OutputMatElastic()
    material.output._configure()
    material.output.writer._configure()

    # Setup integrator
    integrator = ElasticityImplicit()
    integrator.preinitialize(mesh, material)
    return (mesh, integrator)


  def _initialize(self, mesh, integrator):
    """
    Initialize integrator.
    """
    dt = 2.3
    
    from pyre.units.time import s
    integrator.initialize(totalTime=0.0*s, numTimeSteps=1, normalizer=normalizer)
    integrator.timeStep(dt)

    # Setup fields
    from pylith.topology.SolutionFields import SolutionFields
    fields = SolutionFields(mesh)
    fields.add("residual")
    fields.add("disp(t), bc(t+dt)")
    fields.solutionName("disp(t), bc(t+dt)")

    residual = fields.get("residual")
    residual.newSection(residual.VERTICES_FIELD, mesh.coordsys().spaceDim())
    residual.allocate()
    fields.copyLayout("residual")

    residual.zero()
    return


# End of file 
