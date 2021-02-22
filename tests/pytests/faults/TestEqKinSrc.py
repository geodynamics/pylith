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
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file tests/pytests/faults/TestEqKinSrc.py

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
    return


  def test_configure(self):
    """
    Test constructor.
    """
    eqsrc = EqKinSrc()
    eqsrc._configure()
    return


  def test_initialize(self):
    """
    Test initialize().
    """
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    from pythia.pyre.units.time import second

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
    
    from pylith.faults.StepSlipFn import StepSlipFn
    slipfn = StepSlipFn()
    slipfn.inventory.dbSlip = dbFinalSlip
    slipfn.inventory.dbSlipTime = dbSlipTime
    slipfn._configure()

    eqsrc = EqKinSrc()
    eqsrc.inventory.originTime = 5.3*second
    eqsrc.inventory.slipfn = slipfn
    eqsrc._configure()
    eqsrc.preinitialize()
    eqsrc.verifyConfiguration()
    eqsrc.initialize()
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.faults.EqKinSrc import eq_kinematic_src
    fn = eq_kinematic_src()
    return


# End of file 
