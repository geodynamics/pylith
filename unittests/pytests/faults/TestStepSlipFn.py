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

## @file unittests/pytests/faults/TestStepSlipFn.py

## @brief Unit testing of StepSlipFn object.

import unittest

from pylith.faults.StepSlipFn import StepSlipFn

# ----------------------------------------------------------------------
class TestStepSlipFn(unittest.TestCase):
  """
  Unit testing of StepSlipFn object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    slipFn = StepSlipFn()
    slipFn._createCppHandle()
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
    
    slipFn = StepSlipFn()
    slipFn.slip = dbFinalSlip
    slipFn.slipTime = dbSlipTime
    slipFn.preinitialize()
    slipFn.initialize()
    return


# End of file 
