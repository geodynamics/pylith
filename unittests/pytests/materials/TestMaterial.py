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

## @file unittests/pytests/materials/TestMaterial.py

## @brief Unit testing of Material object.

import unittest

# ----------------------------------------------------------------------
class TestMaterial(unittest.TestCase):
  """
  Unit testing of Material object.
  """

  def setUp(self):
    """
    Setup test subject.
    """
    from pylith.materials.ElasticStrain1D import ElasticStrain1D
    self.material = ElasticStrain1D()
    return
    

  def testId(self):
    """
    Test id().
    """
    id = 1234
    self.material.id(id)
    self.assertEqual(id, self.material.id())
    return


  def testLabel(self):
    """
    Test label().
    """
    label = "material abc"
    self.material.label(label)

    # No test of result.
    return


  def testTimeStep(self):
    """
    Test timeStep().
    """
    dt = 0.5
    self.material.timeStep(dt)
    self.assertEqual(dt, self.material.timeStep())
    return


  def testDBProperties(self):
    """
    Test dbProperties().
    """
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    iohandler = SimpleIOAscii()
    iohandler.inventory.filename = "data/matinitialize.spatialdb"
    iohandler._configure()
    db = SimpleDB()
    db.inventory.label = "material properties"
    db.inventory.iohandler = iohandler
    db._configure()

    self.material.dbProperties(db)

    # No test of result.
    return


  def testDBInitialState(self):
    """
    Test dbInitialState().
    """
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    iohandler = SimpleIOAscii()
    iohandler.inventory.filename = "data/matinitialize.spatialdb"
    iohandler._configure()
    db = SimpleDB()
    db.inventory.label = "material properties"
    db.inventory.iohandler = iohandler
    db._configure()

    self.material.dbInitialState(db)

    # No test of result.
    return


  def testNormalizer(self):
    """
    Test normalizer().
    """
    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    normalizer._configure()


    self.material.normalizer(normalizer)

    # No test of result.
    return


  def test_preinitialize(self):
    """
    Test preinitialize().

    WARNING: This is not a rigorous test of initialize() because we
    don't verify the results.
    """
    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATSimplex()
    cell.inventory.shape = "line"
    cell.inventory.order = 1
    cell.inventory.degree = 1
    cell._configure()

    from pylith.feassemble.Quadrature import MeshQuadrature
    quadrature = MeshQuadrature()
    quadrature.inventory.cell = cell
    quadrature.inventory.minJacobian = 1.0e-4
    quadrature._configure()

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    iohandler = SimpleIOAscii()
    iohandler.inventory.filename = "data/matinitialize.spatialdb"
    iohandler._configure()
    db = SimpleDB()
    db.inventory.label = "material properties"
    db.inventory.iohandler = iohandler
    db._configure()

    from pylith.materials.ElasticStrain1D import ElasticStrain1D
    material = ElasticStrain1D()
    material.inventory.quadrature = quadrature
    material.inventory.db = db
    material.inventory.label = "my material"
    material.inventory.id = 54
    material._configure()

    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()
    cs.inventory.spaceDim = 1
    cs._configure()

    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    normalizer._configure()

    from pylith.meshio.MeshIOAscii import MeshIOAscii
    importer = MeshIOAscii()
    importer.inventory.filename = "data/twoelems.mesh"
    importer.inventory.coordsys = cs
    importer._configure()
    mesh = importer.read(normalizer, debug=False, interpolate=False)
    
    material.preinitialize()

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return


# End of file 
