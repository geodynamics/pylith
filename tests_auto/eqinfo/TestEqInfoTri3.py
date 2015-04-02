#!/usr/bin/env python
#
# ----------------------------------------------------------------------
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
# ----------------------------------------------------------------------
#

## @file tests_auto/eqinfo/TestEqInfoTri3.py
##
## @brief Test suite for testing pylith_eqinfo with tri3 fault meshes.

import numpy

from TestEqInfo import TestEqInfo


# Local version of EqInfoApp
from pylith.apps.EqInfoApp import EqInfoApp
class Tri3App(EqInfoApp):
  def __init__(self):
    EqInfoApp.__init__(self, name="tri3")
    return


# Helper function to run pylith_eqinfo.
def run_eqinfo():
  """
  Run pylith_eqinfo.
  """
  if not "done" in dir(run_eqinfo):
    app = Tri3App()
    app.run()
    run_eqinfo.done = True
  return


class TestEqInfoTri3(TestEqInfo):
  """
  Test suite for testing pylith_eqinfo with tri3 meshes.
  """

  def setUp(self):
    """
    Setup for test.
    """
    run_eqinfo()
    return


  def test_stats(self):
    """
    Check fault stats.
    """
    import stats_tri3
    
    timestamp = numpy.array([0.0, 1.0], dtype=numpy.float64)

    oneE = stats_tri3.RuptureStats()
    oneE.timestamp = timestamp
    oneE.ruparea = numpy.array([1.5+2.0, 1.5+2.0], dtype=numpy.float64)
    slip0 = (0.2**2+0.5**2)**0.5
    slip1 = (0.5**2+0.4**2)**0.5
    oneE.potency = numpy.array([slip0*1.5+slip1*2.0, 0.1*1.5+0.2*2.0], dtype=numpy.float64)
    oneE.moment = numpy.array([slip0*1.5*1.0e+10+slip1*2.0*2.0e+10, 
                               0.1*1.5*1.0e+10+0.2*2.0*2.0e+10], dtype=numpy.float64)
    
    twoE = stats_tri3.RuptureStats()
    twoE.timestamp = timestamp
    twoE.ruparea = numpy.array([1.5, 0.0], dtype=numpy.float64)
    twoE.potency = numpy.array([0.1*1.5, 0.0], dtype=numpy.float64)
    twoE.moment = numpy.array([0.1*1.5*1.0e+10, 0.0], dtype=numpy.float64)

    allE = stats_tri3.RuptureStats()
    allE.timestamp = timestamp
    allE.ruparea = oneE.ruparea + twoE.ruparea
    allE.potency = oneE.potency + twoE.potency
    allE.moment = oneE.moment + twoE.moment
    
    self._check(oneE, stats_tri3.one)
    self._check(twoE, stats_tri3.two)
    self._check(allE, stats_tri3.all)
    return


# ----------------------------------------------------------------------
if __name__ == '__main__':
  import unittest
  from TestEqInfoTri3 import TestEqInfoTri3 as Tester

  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(Tester))
  unittest.TextTestRunner(verbosity=2).run(suite)


# End of file 
