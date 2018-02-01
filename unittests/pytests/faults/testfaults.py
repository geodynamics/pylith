#!/usr/bin/env nemesis
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

## @file unittests/faults/testfaults.py

## @brief Python application for testing faults code.

from pylith.tests.UnitTestApp import UnitTestApp

import unittest

class TestApp(UnitTestApp):
  """
  Test application.
  """

  def __init__(self):
    """
    Constructor.
    """
    UnitTestApp.__init__(self)
    return


  def _suite(self):
    """
    Setup the test suite.
    """

    suite = unittest.TestSuite()

    from TestStepSlipFn import TestStepSlipFn
    suite.addTest(unittest.makeSuite(TestStepSlipFn))

    from TestConstRateSlipFn import TestConstRateSlipFn
    suite.addTest(unittest.makeSuite(TestConstRateSlipFn))

    from TestBruneSlipFn import TestBruneSlipFn
    suite.addTest(unittest.makeSuite(TestBruneSlipFn))

    from TestLiuCosSlipFn import TestLiuCosSlipFn
    suite.addTest(unittest.makeSuite(TestLiuCosSlipFn))

    from TestTimeHistorySlipFn import TestTimeHistorySlipFn
    suite.addTest(unittest.makeSuite(TestTimeHistorySlipFn))

    from TestEqKinSrc import TestEqKinSrc
    suite.addTest(unittest.makeSuite(TestEqKinSrc))

    from TestTractPerturbation import TestTractPerturbation
    suite.addTest(unittest.makeSuite(TestTractPerturbation))

    from TestFaultCohesiveKin import TestFaultCohesiveKin
    suite.addTest(unittest.makeSuite(TestFaultCohesiveKin))

    from TestFaultCohesiveDyn import TestFaultCohesiveDyn
    suite.addTest(unittest.makeSuite(TestFaultCohesiveDyn))

    from TestFaultCohesiveImpulses import TestFaultCohesiveImpulses
    suite.addTest(unittest.makeSuite(TestFaultCohesiveImpulses))

    from TestSingleRupture import TestSingleRupture
    suite.addTest(unittest.makeSuite(TestSingleRupture))

    return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = TestApp()
  app.run()


# End of file 
