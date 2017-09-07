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

## @file pylith/tests/FullTestApp.py

## @brief Python application for Python full-scale tests.

import unittest

class FullTestApp(object):
  """
  Test application.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    return


  def main(self):
    """
    Run the application.
    """
    success = unittest.TextTestRunner(verbosity=2).run(self._suite()).wasSuccessful()
    
    if not success:
      import sys
      sys.exit(1)
    return


# End of file 
