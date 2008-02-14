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

## @file unittests/bc/testbc.py

## @brief Python application for testing bc code.

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

    from TestDirichletBoundary import TestDirichletBoundary
    suite.addTest(unittest.makeSuite(TestDirichletBoundary))

    from TestDirichletPoints import TestDirichletPoints
    suite.addTest(unittest.makeSuite(TestDirichletPoints))

    from TestAbsorbingDampers import TestAbsorbingDampers
    suite.addTest(unittest.makeSuite(TestAbsorbingDampers))

    from TestBCSingle import TestBCSingle
    suite.addTest(unittest.makeSuite(TestBCSingle))

    from TestBCTwoSides import TestBCTwoSides
    suite.addTest(unittest.makeSuite(TestBCTwoSides))

    from TestBCFourSides import TestBCFourSides
    suite.addTest(unittest.makeSuite(TestBCFourSides))

    from TestBCSixSides import TestBCSixSides
    suite.addTest(unittest.makeSuite(TestBCSixSides))

    from TestNeumann import TestNeumann
    suite.addTest(unittest.makeSuite(TestNeumann))

    return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = TestApp()
  app.run()


# End of file 
