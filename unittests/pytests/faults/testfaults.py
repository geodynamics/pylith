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
# Copyright (c) 2010-2014 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file unittests/faults/testfaults.py

## @brief Python application for testing faults code.

from pyre.applications.Script import Script

import unittest

class TestApp(Script):
  """
  Test application.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="testapp"):
    """
    Constructor.
    """
    Script.__init__(self, name)
    return


  def main(self):
    """
    Run the application.
    """
    from pylith.utils.PetscManager import PetscManager
    petsc = PetscManager()
    petsc.options = [("malloc_dump", "true")]
    petsc.initialize()

    unittest.TextTestRunner(verbosity=2).run(self._suite())

    petsc.finalize()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

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
