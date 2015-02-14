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

## @file unittests/pytests/faults/TestTimeHistorySlipFn.py

## @brief Unit testing of TimeHistorySlipFn object.

import unittest

from pylith.faults.TimeHistorySlipFn import TimeHistorySlipFn

# ----------------------------------------------------------------------
class TestTimeHistorySlipFn(unittest.TestCase):
  """
  Unit testing of TimeHistorySlipFn object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    slipFn = TimeHistorySlipFn()
    return


  def test_configure(self):
    """
    Test constructor.
    """
    slipFn = TimeHistorySlipFn()
    slipFn._configure()
    return


  def test_initialize(self):
    """
    Test initialize().
    """
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    from spatialdata.spatialdb.TimeHistory import TimeHistory

    ioFinalSlip = SimpleIOAscii()
    ioFinalSlip.inventory.filename = "finalslip.spatialdb"
    ioFinalSlip._configure()
    dbFinalSlip = SimpleDB()
    dbFinalSlip.inventory.iohandler = ioFinalSlip
    dbFinalSlip.inventory.label = "slip amplitude"
    dbFinalSlip._configure()
    
    ioSlipTime = SimpleIOAscii()
    ioSlipTime.inventory.filename = "sliptime.spatialdb"
    ioSlipTime._configure()
    dbSlipTime = SimpleDB()
    dbSlipTime.inventory.iohandler = ioSlipTime
    dbSlipTime.inventory.label = "slip time"
    dbSlipTime._configure()
    
    dbTimeHistory = TimeHistory()
    dbTimeHistory.inventory.filename = "slipfn.timedb"
    dbTimeHistory.inventory.label = "time history"
    dbTimeHistory._configure()
    
    slipFn = TimeHistorySlipFn()
    slipFn.inventory.dbslip = dbFinalSlip
    slipFn.inventory.dbSlipTime = dbSlipTime
    slipFn.inventory.dbTimeHistory = dbTimeHistory
    slipFn._configure()
    slipFn.preinitialize()
    slipFn.verifyConfiguration()
    slipFn.initialize()
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.faults.TimeHistorySlipFn import slip_time_fn
    fn = slip_time_fn()
    return


# End of file 
