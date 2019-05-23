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
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================

from pylith.tests.FullTestApp import FullTestApp

import unittest

class TestApp(FullTestApp):
  """
  Test application.
  """

  def __init__(self):
    """
    Constructor.
    """
    FullTestApp.__init__(self)
    return


  def _suite(self):
    """
    Create test suite.
    """
    suite = unittest.TestSuite()

    from TestAxialDisp import TestAxialDisp
    suite.addTest(unittest.makeSuite(TestAxialDisp))

    #from TestShearTraction import TestShearTraction
    #suite.addTest(unittest.makeSuite(TestShearTraction))

    return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
  TestApp().main()

  
# End of file 
