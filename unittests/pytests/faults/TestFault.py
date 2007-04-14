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

## @file unittests/pytests/faults/TestFault.py

## @brief Unit testing of Fault object.

import unittest

from pylith.faults.Fault import Fault

# ----------------------------------------------------------------------
class TestFault(unittest.TestCase):
  """
  Unit testing of Fault object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    fault = Fault()
    self.assertEqual(None, fault.cppHandle)
    return


  def test_initialize(self):
    """
    Test initialize().
    """
    raise NotImplementedError("Need to implement unit test.")
    return
  

# End of file 
