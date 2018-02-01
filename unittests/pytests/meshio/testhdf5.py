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
#

## @file unittests/meshio/testhdf5.py

## @brief Python application for testing HDF5 related meshio code.

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

    from TestDataWriterHDF5 import TestDataWriterHDF5
    suite.addTest(unittest.makeSuite(TestDataWriterHDF5))

    from TestDataWriterHDF5Ext import TestDataWriterHDF5Ext
    suite.addTest(unittest.makeSuite(TestDataWriterHDF5Ext))

    from TestXdmf import TestXdmf
    suite.addTest(unittest.makeSuite(TestXdmf))

    return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = TestApp()
  app.run()


# End of file 
