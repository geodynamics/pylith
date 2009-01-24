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

## @file unittests/meshio/testmeshio.py

## @brief Python application for testing meshio code.

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

    from TestMeshIOAscii import TestMeshIOAscii
    suite.addTest(unittest.makeSuite(TestMeshIOAscii))

    #from TestMeshIOLagrit import TestMeshIOLagrit
    #suite.addTest(unittest.makeSuite(TestMeshIOLagrit))

    #from TestVertexFilterVecNorm import TestVertexFilterVecNorm
    #suite.addTest(unittest.makeSuite(TestVertexFilterVecNorm))

    #from TestCellFilterAvg import TestCellFilterAvg
    #suite.addTest(unittest.makeSuite(TestCellFilterAvg))

    #from TestOutputManager import TestOutputManager
    #suite.addTest(unittest.makeSuite(TestOutputManager))

    #from TestOutputSolnSubset import TestOutputSolnSubset
    #suite.addTest(unittest.makeSuite(TestOutputSolnSubset))

    return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = TestApp()
  app.run()


# End of file 
