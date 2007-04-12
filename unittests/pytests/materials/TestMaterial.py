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

  def test_initialize(self):
    """
    Test initialize().

    WARNING: This is not a rigorous test of initialize() because we
    don't verify the results.
    """
    from pylith.materials.ElasticIsotropic3D import ElasticIsotropic3D
    material = ElasticIsotropic3D()

    from pylith.feassemble.quadrature.Quadrature1D import Quadrature1D
    quadrature = Quadrature1D()
    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATSimplex()
    cell.shape = "line"
    cell.order = 1
    cell.degree = 1
    quadrature.cell = cell
    quadrature.minJacobian = 1.0e-4
    quadrature.initialize()
    material.quadrature = quadrature

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    iohandler = SimpleIOAscii()
    iohandler.filename = "data/matinitialize.spatialdb"
    db = SimpleDB()
    db.label = "material properties"
    db.iohandler = iohandler
    material.db = db
    material.label = "my material"
    material.id = 54

    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()
    cs.spaceDim = 1

    from pylith.meshio.MeshIOAscii import MeshIOAscii
    importer = MeshIOAscii()
    importer.filename = "data/twoelems.mesh"
    importer.coordsys = cs
    mesh = importer.read(debug=False, interpolate=False)
    
    material.initialize(mesh)

    # We should really add something here to check to make sure things
    # actually initialized correctly    
    return


# End of file 
