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

## @file unittests/pytests/materials/TestElasticMaterial.py

## @brief Unit testing of Material object.

import unittest

# ----------------------------------------------------------------------
class TestElasticMaterial(unittest.TestCase):
  """
  Unit testing of ElasticMaterial object.
  """

  def setUp(self):
    """
    Setup test subject.
    """
    from pylith.materials.ElasticStrain1D import ElasticStrain1D
    self.material = ElasticStrain1D()
    return
    

  def testDBInitialStress(self):
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    iohandler = SimpleIOAscii()
    iohandler.inventory.filename = "data/matinitialize.spatialdb"
    iohandler._configure()
    db = SimpleDB()
    db.inventory.label = "material properties"
    db.inventory.iohandler = iohandler
    db._configure()

    self.material.dbInitialStress(db)

    # No test of result.
    return


  def testDBInitialStrain(self):
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    iohandler = SimpleIOAscii()
    iohandler.inventory.filename = "data/matinitialize.spatialdb"
    iohandler._configure()
    db = SimpleDB()
    db.inventory.label = "material properties"
    db.inventory.iohandler = iohandler
    db._configure()

    self.material.dbInitialStrain(db)

    # No test of result.
    return


# End of file 
