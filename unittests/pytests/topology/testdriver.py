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

## @file unittests/topology/testdriver.py

## @brief Python application for testing topology code.

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

    from TestMesh import TestMesh
    suite.addTest(unittest.makeSuite(TestMesh))

    from TestSubMesh import TestSubMesh
    suite.addTest(unittest.makeSuite(TestSubMesh))

    from TestFieldBase import TestFieldBase
    suite.addTest(unittest.makeSuite(TestFieldBase))

    from TestMeshField import TestMeshField
    suite.addTest(unittest.makeSuite(TestMeshField))

    from TestMeshFields import TestMeshFields
    suite.addTest(unittest.makeSuite(TestMeshFields))

    from TestSolutionFields import TestSolutionFields
    suite.addTest(unittest.makeSuite(TestSolutionFields))

    from TestMeshGenerator import TestMeshGenerator
    suite.addTest(unittest.makeSuite(TestMeshGenerator))

    from TestMeshImporter import TestMeshImporter
    suite.addTest(unittest.makeSuite(TestMeshImporter))

    return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = TestApp()
  app.run()


# End of file 
