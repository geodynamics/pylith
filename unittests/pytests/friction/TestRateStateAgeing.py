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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file unittests/pytests/friction/TestRateStateAgeing.py

## @brief Unit testing of RateStateAgeing object.

import unittest

from pylith.friction.RateStateAgeing import RateStateAgeing

# ----------------------------------------------------------------------
class TestRateStateAgeing(unittest.TestCase):
  """
  Unit testing of RateStateAgeing object.
  """

  def setUp(self):
    """
    Setup test subject.
    """
    self.friction = RateStateAgeing()
    return
  

  def test_constructor(self):
    """
    Test constructor.
    """
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.friction.RateStateAgeing import friction_model
    m = friction_model()
    return


# End of file 
