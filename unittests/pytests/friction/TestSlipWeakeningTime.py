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

## @file unittests/pytests/friction/TestSlipWeakeningTime.py

## @brief Unit testing of SlipWeakeningTime object.

import unittest

from pylith.friction.SlipWeakeningTime import SlipWeakeningTime

# ----------------------------------------------------------------------
class TestSlipWeakeningTime(unittest.TestCase):
  """
  Unit testing of SlipWeakeningTime object.
  """

  def setUp(self):
    """
    Setup test subject.
    """
    self.friction = SlipWeakeningTime()
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
    from pylith.friction.SlipWeakeningTime import friction_model
    m = friction_model()
    return


# End of file 
