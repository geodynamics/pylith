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
# Copyright (c) 2010-2011 University of California, Davis
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

    from TestMeshIOLagrit import TestMeshIOLagrit
    suite.addTest(unittest.makeSuite(TestMeshIOLagrit))

    from TestVertexFilterVecNorm import TestVertexFilterVecNormMesh
    suite.addTest(unittest.makeSuite(TestVertexFilterVecNormMesh))

    from TestVertexFilterVecNorm import TestVertexFilterVecNormSubMesh
    suite.addTest(unittest.makeSuite(TestVertexFilterVecNormSubMesh))

    from TestCellFilterAvg import TestCellFilterAvgMesh
    suite.addTest(unittest.makeSuite(TestCellFilterAvgMesh))

    from TestCellFilterAvg import TestCellFilterAvgSubMesh
    suite.addTest(unittest.makeSuite(TestCellFilterAvgSubMesh))

    from TestDataWriterVTK import TestDataWriterVTKMesh
    suite.addTest(unittest.makeSuite(TestDataWriterVTKMesh))

    from TestDataWriterVTK import TestDataWriterVTKSubMesh
    suite.addTest(unittest.makeSuite(TestDataWriterVTKSubMesh))

    from TestDataWriterVTK import TestDataWriterVTKSubSubMesh
    suite.addTest(unittest.makeSuite(TestDataWriterVTKSubSubMesh))

    from TestOutputManagerMesh import TestOutputManagerMesh
    suite.addTest(unittest.makeSuite(TestOutputManagerMesh))

    from TestOutputManagerSubMesh import TestOutputManagerSubMesh
    suite.addTest(unittest.makeSuite(TestOutputManagerSubMesh))

    from TestOutputSolnSubset import TestOutputSolnSubset
    suite.addTest(unittest.makeSuite(TestOutputSolnSubset))

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
