#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2012 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file tests/2d/tri3/TestDislocation.py
##
## @brief Test suite for testing pylith with shear dislocation for 2-D
## box.

import unittest

# Local application
from pylith.apps.PyLithApp import PyLithApp
class DislocationApp(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="dislocation")
    return


# Helper function to run PyLith
def run_pylith():
  """
  Run pylith.
  """
  if not "done" in dir(run_pylith):
    app = DislocationApp()
    app.run()
    run_pylith.done = True
  return


class TestDislocation(unittest.TestCase):
  """
  Test suite for testing pylith with shear dislocation for 2-D box.
  """

  def setUp(self):
    """
    Setup for test.
    """
    run_pylith()
    return


  def test_disp(self):
    """
    Check displacement field.
    """
    return


# ----------------------------------------------------------------------
if __name__ == '__main__':
  import unittest
  from TestDislocation import TestDislocation as Tester

  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(Tester))
  unittest.TextTestRunner(verbosity=2).run(suite)


# End of file 
