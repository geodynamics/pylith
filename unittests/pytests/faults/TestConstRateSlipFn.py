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
    slipFn._createCppHandle()
    self.failIfEqual(None, slipFn.cppHandle)
    return


  def test_initialize(self):
    """
    Test initialize().
    """
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii

    ioSlipRate = SimpleIOAscii()
    ioSlipRate.filename = "sliprate.spatialdb"
    dbSlipRate = SimpleDB()
    dbSlipRate.iohandler = ioSlipRate
    dbSlipRate.label = "slip rate"
    
    ioSlipTime = SimpleIOAscii()
    ioSlipTime.filename = "sliptime.spatialdb"
    dbSlipTime = SimpleDB()
    dbSlipTime.iohandler = ioSlipTime
    dbSlipTime.label = "slip time"
    
    slipFn = ConstRateSlipFn()
    slipFn.slipRate = dbSlipRate
    slipFn.slipTime = dbSlipTime
    slipFn.preinitialize()
    slipFn.initialize()
    return


# End of file 
