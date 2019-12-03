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

## @file tests/pytests/friction/TestTimeWeakening.py

## @brief Unit testing of TimeWeakening object.

import unittest

from pylith.friction.TimeWeakening import TimeWeakening

# ----------------------------------------------------------------------
class TestTimeWeakening(unittest.TestCase):
  """
  Unit testing of TimeWeakening object.
  """

  def setUp(self):
    """
    Setup test subject.
    """
    self.friction = TimeWeakening()
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
    from pylith.friction.TimeWeakening import friction_model
    m = friction_model()
    return


# End of file 
