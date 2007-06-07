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

    from TestBruneSlipFn import TestBruneSlipFn
    suite.addTest(unittest.makeSuite(TestBruneSlipFn))

    from TestFault import TestFault
    suite.addTest(unittest.makeSuite(TestFault))

    from TestFaultCohesive import TestFaultCohesive
    suite.addTest(unittest.makeSuite(TestFaultCohesive))

    from TestFaultCohesiveKin import TestFaultCohesiveKin
    suite.addTest(unittest.makeSuite(TestFaultCohesiveKin))

    return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = TestApp()
  app.run()


# End of file 
