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

## @file unittests/friction/testfriction.py

## @brief Python application for testing fault constitutive models.

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

    from TestStaticFriction import TestStaticFriction
    suite.addTest(unittest.makeSuite(TestStaticFriction))

    from TestSlipWeakening import TestSlipWeakening
    suite.addTest(unittest.makeSuite(TestSlipWeakening))

    from TestSlipWeakeningTime import TestSlipWeakeningTime
    suite.addTest(unittest.makeSuite(TestSlipWeakeningTime))

    from TestSlipWeakeningTimeStable import TestSlipWeakeningTimeStable
    suite.addTest(unittest.makeSuite(TestSlipWeakeningTimeStable))

    from TestRateStateAgeing import TestRateStateAgeing
    suite.addTest(unittest.makeSuite(TestRateStateAgeing))

    from TestTimeWeakening import TestTimeWeakening
    suite.addTest(unittest.makeSuite(TestTimeWeakening))

    from TestFrictionModel import TestFrictionModel
    suite.addTest(unittest.makeSuite(TestFrictionModel))

    return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = TestApp()
  app.run()


# End of file 
