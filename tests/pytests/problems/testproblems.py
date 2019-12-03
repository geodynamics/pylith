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

## @file tests/problems/testproblems.py

## @brief Python application for testing problems code.

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

    from TestSolutionSubfields import TestSolutionSubfields
    suite.addTest(unittest.makeSuite(TestSolutionSubfields))

    from TestSolution import TestSolution
    suite.addTest(unittest.makeSuite(TestSolution))

    from TestTimeStep import TestTimeStep
    suite.addTest(unittest.makeSuite(TestTimeStep))

    from TestTimeStepUniform import TestTimeStepUniform
    suite.addTest(unittest.makeSuite(TestTimeStepUniform))

    from TestTimeStepUser import TestTimeStepUser
    suite.addTest(unittest.makeSuite(TestTimeStepUser))

    from TestTimeStepAdapt import TestTimeStepAdapt
    suite.addTest(unittest.makeSuite(TestTimeStepAdapt))

    from TestProgressMonitor import TestProgressMonitor
    suite.addTest(unittest.makeSuite(TestProgressMonitor))

    from TestProgressMonitorTime import TestProgressMonitorTime
    suite.addTest(unittest.makeSuite(TestProgressMonitorTime))

    from TestProgressMonitorStep import TestProgressMonitorStep
    suite.addTest(unittest.makeSuite(TestProgressMonitorStep))

    return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = TestApp()
  app.run()


# End of file 
