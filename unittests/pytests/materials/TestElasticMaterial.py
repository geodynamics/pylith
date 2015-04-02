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
    from pylith.materials.ElasticPlaneStrain import ElasticPlaneStrain
    self.material = ElasticPlaneStrain()
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

    material = self.material
    material.inventory.label = "Elastic material"
    material.inventory.dbInitialStress = db
    material.inventory.useInitialStress = True
    material._configure()

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

    material = self.material
    material.inventory.label = "Elastic material"
    material.inventory.dbInitialStrain = db
    material.inventory.useInitialStrain = True
    material._configure()

    # No test of result.
    return


# End of file 
