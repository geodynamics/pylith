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

## @file unittests/pytests/faults/TestEqKinSrc.py

## @brief Unit testing of EqKinSrc object.

import unittest

from pylith.faults.EqKinSrc import EqKinSrc

# ----------------------------------------------------------------------
class TestEqKinSrc(unittest.TestCase):
  """
  Unit testing of EqKinSrc object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    eqsrc = EqKinSrc()
    eqsrc._createCppHandle()
    self.failIfEqual(None, eqsrc.cppHandle)
    return


  def test_initialize(self):
    """
    Test initialize().
    """
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    from pyre.units.time import second

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
    
    from pylith.faults.BruneSlipFn import BruneSlipFn
    slipfn = BruneSlipFn()
    slipfn.slip = dbFinalSlip
    slipfn.slipTime = dbSlipTime
    slipfn.slipRate = dbPeakRate

    eqsrc = EqKinSrc()
    eqsrc.originTime = 5.3*second
    eqsrc.slipfn = slipfn
    eqsrc.preinitialize()
    eqsrc.initialize()
    return


# End of file 
