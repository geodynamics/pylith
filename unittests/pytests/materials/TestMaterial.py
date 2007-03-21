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
    from pylith.utils.PetscManager import PetscManager
    petsc = PetscManager()
    petsc.initialize()

    from pylith.materials.ElasticIsotropic3D import ElasticIsotropic3D
    material = ElasticIsotropic3D()

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    iohandler = SimpleIOAscii()
    iohandler.filename = "data/matinitialize.spatialdb"
    db = SimpleDB()
    db.iohandler = iohandler
    material.db = db
    material.label = "my material"
    material.id = 54

    from pylith.meshio.MeshIOAscii import MeshIOAscii
    importer = MeshIOAscii()
    importer.filename = "data/twoelems.mesh"
    mesh = importer.read()
    
    material.initialize(mesh)

    # We should really add something here to check to make sure things
    # actually initialized correctly    

    petsc.finalize()
    return


# End of file 
