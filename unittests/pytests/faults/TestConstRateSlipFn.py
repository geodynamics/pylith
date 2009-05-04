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

## @file unittests/pytests/faults/TestConstRateSlipFn.py

## @brief Unit testing of ConstRateSlipFn object.

import unittest

from pylith.faults.ConstRateSlipFn import ConstRateSlipFn

# ----------------------------------------------------------------------
class TestConstRateSlipFn(unittest.TestCase):
  """
  Unit testing of ConstRateSlipFn object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    slipFn = ConstRateSlipFn()
    return


  def test_configure(self):
    """
    Test _configure().
    """
    slipFn = ConstRateSlipFn()
    slipFn._configure()
    return


  def test_initialize(self):
    """
    Test initialize().
    """
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii

    ioSlipRate = SimpleIOAscii()
    ioSlipRate.inventory.filename = "sliprate.spatialdb"
    ioSlipRate._configure()
    dbSlipRate = SimpleDB()
    dbSlipRate.inventory.iohandler = ioSlipRate
    dbSlipRate.inventory.label = "slip rate"
    dbSlipRate._configure()
    
    ioSlipTime = SimpleIOAscii()
    ioSlipTime.inventory.filename = "sliptime.spatialdb"
    ioSlipTime._configure()
    dbSlipTime = SimpleDB()
    dbSlipTime.inventory.iohandler = ioSlipTime
    dbSlipTime.inventory.label = "slip time"
    dbSlipTime._configure()
    
    slipFn = ConstRateSlipFn()
    slipFn.inventory.dbSlipRate = dbSlipRate
    slipFn.inventory.dbSlipTime = dbSlipTime
    slipFn._configure()
    slipFn.preinitialize()
    slipFn.verifyConfiguration()
    slipFn.initialize()
    return


# End of file 
