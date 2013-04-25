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
# Copyright (c) 2010-2013 University of California, Davis
#
# See COPYING for license information.
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

    from TestCellGeometry import TestCellGeometry
    suite.addTest(unittest.makeSuite(TestCellGeometry))

    from TestFIATSimplex import TestFIATSimplex
    suite.addTest(unittest.makeSuite(TestFIATSimplex))

    from TestFIATLagrange import TestFIATLagrange
    suite.addTest(unittest.makeSuite(TestFIATLagrange))

    from TestMeshQuadrature import TestMeshQuadrature
    suite.addTest(unittest.makeSuite(TestMeshQuadrature))

    from TestSubMeshQuadrature import TestSubMeshQuadrature
    suite.addTest(unittest.makeSuite(TestSubMeshQuadrature))

    from TestElasticityImplicit import TestElasticityImplicit
    suite.addTest(unittest.makeSuite(TestElasticityImplicit))

    from TestElasticityExplicit import TestElasticityExplicit
    suite.addTest(unittest.makeSuite(TestElasticityExplicit))

    from TestElasticityImplicitLgDeform import TestElasticityImplicitLgDeform
    suite.addTest(unittest.makeSuite(TestElasticityImplicitLgDeform))

    from TestElasticityExplicitLgDeform import TestElasticityExplicitLgDeform
    suite.addTest(unittest.makeSuite(TestElasticityExplicitLgDeform))

    return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = TestApp()
  app.run()


# End of file 
