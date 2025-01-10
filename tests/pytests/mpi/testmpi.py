#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file tests/bc/testmpi.py

# @brief Python application for testing mpi code.

from pylith.tests.UnitTestApp import UnitTestApp

import unittest

class TestApp(UnitTestApp):
  """
  Test application.
  """

  def __init__(self):
    """
    Constructor.
    """
    UnitTestApp.__init__(self)
    return


  def _suite(self):
    """
    Setup the test suite.
    """

    suite = unittest.TestSuite()

    from TestCommunicator import TestCommunicator
    suite.addTest(unittest.makeSuite(TestCommunicator))

    from TestReduce import TestReduce
    suite.addTest(unittest.makeSuite(TestReduce))

    return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = TestApp()
  app.run()


# End of file 
