#!/usr/bin/env python
#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# {LicenseText}
#
# ======================================================================
#

## @file unittests/pytests/friction/TestStaticFriction.py

## @brief Unit testing of StaticFriction object.

import unittest

from pylith.friction.StaticFriction import StaticFriction

# ----------------------------------------------------------------------
class TestStaticFriction(unittest.TestCase):
  """
  Unit testing of StaticFriction object.
  """

  def setUp(self):
    """
    Setup test subject.
    """
    self.friction = StaticFriction()
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
    from pylith.friction.StaticFriction import friction_model
    m = friction_model()
    return


# End of file 
