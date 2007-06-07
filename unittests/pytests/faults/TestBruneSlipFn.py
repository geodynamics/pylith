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

## @file unittests/pytests/faults/TestBruneSlipFn.py

## @brief Unit testing of BruneSlipFn object.

import unittest

from pylith.faults.BruneSlipFn import BruneSlipFn

# ----------------------------------------------------------------------
class TestBruneSlipFn(unittest.TestCase):
  """
  Unit testing of BruneSlipFn object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    slipFn = BruneSlipFn()
    self.failIfEqual(None, slipFn.cppHandle)
    return


  def test_initialize(self):
    """
    Test initialize().
    """
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii

    ioFinalSlip = SimpleIOAscii()
    ioFinalSlip.filename = "finalslip.spatialdb"
    dbFinalSlip = SimpleDB()
    dbFinalSlip.iohandler = ioFinalSlip
    dbFinalSlip.label = "final slip"
    
    ioSlipTime = SimpleIOAscii()
    ioSlipTime.filename = "sliptime.spatialdb"
    dbSlipTime = SimpleDB()
    dbSlipTime.iohandler = ioSlipTime
    dbSlipTime.label = "slip time"
    
    ioPeakRate = SimpleIOAscii()
    ioPeakRate.filename = "peakrate.spatialdb"
    dbPeakRate = SimpleDB()
    dbPeakRate.iohandler = ioPeakRate
    dbPeakRate.label = "peak rate"
    
    slipFn = BruneSlipFn()
    slipFn.slip = dbFinalSlip
    slipFn.slipTime = dbSlipTime
    slipFn.slipRate = dbPeakRate
    slipFn.initialize()
    return


# End of file 
