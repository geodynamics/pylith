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
    return


  def test_configure(self):
    """
    Test constructor.
    """
    slipFn = StepSlipFn()
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
    
    slipFn = StepSlipFn()
    slipFn.inventory.dbSlip = dbFinalSlip
    slipFn.inventory.dbSlipTime = dbSlipTime
    slipFn._configure()
    slipFn.preinitialize()
    slipFn.verifyConfiguration()
    slipFn.initialize()
    return


# End of file 
