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

## @file unittests/pytests/utils/TestPylith.py

## @brief Unit testing of Pylith module.

import unittest


# ----------------------------------------------------------------------
class TestPylith(unittest.TestCase):
  """
  Unit testing of Pylith object.
  """
  

  def test_sizeofPylithScalar(self):
    """
    Test constructor.
    """
    from pylith.utils.utils import sizeofPylithScalar
    size = sizeofPylithScalar()
    self.failUnless(4 == size or 8 == size)
    return


# End of file 
