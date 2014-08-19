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
# Copyright (c) 2010-2014 University of California, Davis
#
# See COPYING for license information.
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

    from TestJacobian import TestJacobian
    suite.addTest(unittest.makeSuite(TestJacobian))

    from TestMeshGenerator import TestMeshGenerator
    suite.addTest(unittest.makeSuite(TestMeshGenerator))

    from TestMeshImporter import TestMeshImporter
    suite.addTest(unittest.makeSuite(TestMeshImporter))

    #print "WARNING: TestRefineUniform NOT IMPLEMENTED"
    from TestRefineUniform import TestRefineUniform
    suite.addTest(unittest.makeSuite(TestRefineUniform))

    return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = TestApp()
  app.run()


# End of file 
