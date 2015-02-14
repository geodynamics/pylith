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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
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

    from TestMeshIOAscii import TestMeshIOAscii
    suite.addTest(unittest.makeSuite(TestMeshIOAscii))

    from TestMeshIOLagrit import TestMeshIOLagrit
    suite.addTest(unittest.makeSuite(TestMeshIOLagrit))

    from TestVertexFilterVecNorm import TestVertexFilterVecNorm
    suite.addTest(unittest.makeSuite(TestVertexFilterVecNorm))

    from TestCellFilterAvg import TestCellFilterAvg
    suite.addTest(unittest.makeSuite(TestCellFilterAvg))

    from TestDataWriterVTK import TestDataWriterVTK
    suite.addTest(unittest.makeSuite(TestDataWriterVTK))

    from TestOutputManagerMesh import TestOutputManagerMesh
    suite.addTest(unittest.makeSuite(TestOutputManagerMesh))

    from TestOutputManagerSubMesh import TestOutputManagerSubMesh
    suite.addTest(unittest.makeSuite(TestOutputManagerSubMesh))

    from TestOutputSolnSubset import TestOutputSolnSubset
    suite.addTest(unittest.makeSuite(TestOutputSolnSubset))

    from TestOutputSolnPoints import TestOutputSolnPoints
    suite.addTest(unittest.makeSuite(TestOutputSolnPoints))

    from TestSingleOutput import TestSingleOutput
    suite.addTest(unittest.makeSuite(TestSingleOutput))

    #TestOutputNeumann

    #TestOutputFaultKin

    #TestOutputDirichlet

    return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = TestApp()
  app.run()


# End of file 
