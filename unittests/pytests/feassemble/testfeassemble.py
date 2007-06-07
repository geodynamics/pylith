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

## @file unittests/feassemble/testfeassemble.py

## @brief Python application for testing feassemble code.

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
    manager = PetscManager()
    manager.initialize()
    unittest.TextTestRunner(verbosity=2).run(self._suite())
    manager.finalize()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _suite(self):
    """
    Setup the test suite.
    """

    suite = unittest.TestSuite()

    from TestCellGeometry import TestCellGeometry
    suite.addTest(unittest.makeSuite(TestCellGeometry))

    from TestFIATSimplex import TestFIATSimplex
    suite.addTest(unittest.makeSuite(TestFIATSimplex))

    from TestFIATLagrange import TestFIATLagrange
    suite.addTest(unittest.makeSuite(TestFIATLagrange))

    from TestQuadrature import TestQuadrature
    suite.addTest(unittest.makeSuite(TestQuadrature))

    from TestIntegrator import TestIntegrator
    suite.addTest(unittest.makeSuite(TestIntegrator))

    from TestElasticityExplicit import TestElasticityExplicit
    suite.addTest(unittest.makeSuite(TestElasticityExplicit))

    from TestElasticityImplicit import TestElasticityImplicit
    suite.addTest(unittest.makeSuite(TestElasticityImplicit))

    return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = TestApp()
  app.run()


# End of file 
