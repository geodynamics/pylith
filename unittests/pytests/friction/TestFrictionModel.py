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

## @file unittests/pytests/friction/TestFrictionModel.py

## @brief Unit testing of Material object.

import unittest

# ----------------------------------------------------------------------
class TestFrictionModel(unittest.TestCase):
  """
  Unit testing of Material object.
  """

  def setUp(self):
    """
    Setup test subject.
    """
    from pylith.friction.StaticFriction import StaticFriction
    self.friction = StaticFriction()
    return
    

  def testLabel(self):
    """
    Test label().
    """
    label = "friction abc"
    self.friction.label(label)
    self.assertEqual(label, self.friction.label())
    return


  def testTimeStep(self):
    """
    Test timeStep().
    """
    dt = 0.5
    self.friction.timeStep(dt)
    self.assertEqual(dt, self.friction.timeStep())
    return


  def testDBProperties(self):
    """
    Test dbProperties().
    """
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    iohandler = SimpleIOAscii()
    iohandler.inventory.filename = "data/staticfriction.spatialdb"
    iohandler._configure()
    db = SimpleDB()
    db.inventory.label = "friction properties"
    db.inventory.iohandler = iohandler
    db._configure()

    self.friction.dbProperties(db)

    # No test of result.
    return


  def testDBInitialState(self):
    """
    Test dbInitialState().
    """
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    iohandler = SimpleIOAscii()
    iohandler.inventory.filename = "data/staticfriction.spatialdb"
    iohandler._configure()
    db = SimpleDB()
    db.inventory.label = "friction properties"
    db.inventory.iohandler = iohandler
    db._configure()

    self.friction.dbInitialState(db)

    # No test of result.
    return


  def testNormalizer(self):
    """
    Test normalizer().
    """
    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    normalizer._configure()


    self.friction.normalizer(normalizer)

    # No test of result.
    return


# End of file 
