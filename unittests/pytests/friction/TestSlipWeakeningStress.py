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

## @file unittests/pytests/friction/TestSlipWeakeningStress.py

## @brief Unit testing of SlipWeakeningStress object.

import unittest

from pylith.friction.SlipWeakeningStress import SlipWeakeningStress

# ----------------------------------------------------------------------
class TestSlipWeakeningStress(unittest.TestCase):
  """
  Unit testing of SlipWeakeningStress object.
  """

  def setUp(self):
    """
    Setup test subject.
    """
    self.friction = SlipWeakeningStress()
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
    from pylith.friction.SlipWeakeningStress import friction_model
    m = friction_model()
    return


# End of file 
