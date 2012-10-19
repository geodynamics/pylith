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
# Copyright (c) 2010-2012 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file tests_auto/eqinfo/TestEqInfoQuad4.py
##
## @brief Test suite for testing pylith_eqinfo with quad4 fault meshes.

import numpy

from TestEqInfo import TestEqInfo


# Local version of EqInfoApp
from pylith.apps.EqInfoApp import EqInfoApp
class Quad4App(EqInfoApp):
  def __init__(self):
    EqInfoApp.__init__(self, name="quad4")
    return


# Helper function to run pylith_eqinfo.
def run_eqinfo():
  """
  Run pylith_eqinfo.
  """
  if not "done" in dir(run_eqinfo):
    app = Quad4App()
    app.run()
    run_eqinfo.done = True
  return


class TestEqInfoQuad4(TestEqInfo):
  """
  Test suite for testing pylith_eqinfo with quad4 meshes.
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
    import stats_quad4
    
    timestamp = numpy.array([5.0], dtype=numpy.float64)

    area0 = 1.5*1.75
    area1 = 1.35*1.5
    slip0 = (1.0**2+1.2**2)**0.5
    slip1 = (1.4**2+1.6**2)**0.5
    oneE = stats_quad4.RuptureStats()
    oneE.timestamp = timestamp
    oneE.ruparea = numpy.array([area0 + area1], dtype=numpy.float64)
    oneE.potency = numpy.array([slip0*area0 + slip1*area1], dtype=numpy.float64)
    oneE.moment = numpy.array([slip0*area0*1.0e+10 + area1*slip1*2.0e+10], dtype=numpy.float64)
    
    area0 = 1.5*1.75
    area1 = 1.35*1.5
    slip0 = (1.0**2+1.2**2)**0.5
    slip1 = (1.4**2+1.6**2)**0.5
    twoE = stats_quad4.RuptureStats()
    twoE.timestamp = timestamp
    twoE.ruparea = numpy.array([area0 + area1], dtype=numpy.float64)
    twoE.potency = numpy.array([slip0*area0 + slip1*area1], dtype=numpy.float64)
    twoE.moment = numpy.array([slip0*area0*1.0e+10 + area1*slip1*2.0e+10], dtype=numpy.float64)

    allE = stats_quad4.RuptureStats()
    allE.timestamp = timestamp
    allE.ruparea = oneE.ruparea + twoE.ruparea
    allE.potency = oneE.potency + twoE.potency
    allE.moment = oneE.moment + twoE.moment
    
    self._check(oneE, stats_quad4.one)
    self._check(twoE, stats_quad4.two)
    self._check(allE, stats_quad4.all)
    return


# End of file 
