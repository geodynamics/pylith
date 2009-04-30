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

## @file unittests/pytests/faults/TestLiuCosSlipFn.py

## @brief Unit testing of LiuCosSlipFn object.

import unittest

from pylith.faults.LiuCosSlipFn import LiuCosSlipFn

# ----------------------------------------------------------------------
class TestLiuCosSlipFn(unittest.TestCase):
  """
  Unit testing of LiuCosSlipFn object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    slipFn = LiuCosSlipFn()
    return


  def test_configure(self):
    """
    Test constructor.
    """
    slipFn = LiuCosSlipFn()
    slipFn._configure()
    return


  def test_initialize(self):
    """
    Test initialize().
    """
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii

    ioFinalSlip = SimpleIOAscii()
    ioFinalSlip.inventory.filename = "finalslip.spatialdb"
    ioFinalSlip._configure()
    dbFinalSlip = SimpleDB()
    dbFinalSlip.inventory.iohandler = ioFinalSlip
    dbFinalSlip.inventory.label = "final slip"
    dbFinalSlip._configure()
    
    ioSlipTime = SimpleIOAscii()
    ioSlipTime.inventory.filename = "sliptime.spatialdb"
    ioSlipTime._configure()
    dbSlipTime = SimpleDB()
    dbSlipTime.inventory.iohandler = ioSlipTime
    dbSlipTime.inventory.label = "slip time"
    dbSlipTime._configure()
    
    ioRiseTime = SimpleIOAscii()
    ioRiseTime.inventory.filename = "risetime.spatialdb"
    ioRiseTime._configure()
    dbRiseTime = SimpleDB()
    dbRiseTime.inventory.iohandler = ioRiseTime
    dbRiseTime.inventory.label = "rise time"
    dbRiseTime._configure()
    
    slipFn = LiuCosSlipFn()
    slipFn.inventory.dbslip = dbFinalSlip
    slipFn.inventory.dbSlipTime = dbSlipTime
    slipFn.inventory.dbRiseTime = dbRiseTime
    slipFn._configure()
    slipFn.preinitialize()
    slipFn.verifyConfiguration()
    slipFn.initialize()
    return


# End of file 
