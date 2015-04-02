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


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.faults.ConstRateSlipFn import slip_time_fn
    fn = slip_time_fn()
    return


# End of file 
