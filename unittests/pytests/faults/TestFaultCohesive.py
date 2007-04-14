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

## @file unittests/pytests/faults/TestFaultCohesive.py

## @brief Unit testing of Fault object.

import unittest

from pylith.faults.FaultCohesive import FaultCohesive

# ----------------------------------------------------------------------
class TestFaultCohesive(unittest.TestCase):
  """
  Unit testing of Fault object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    fault = FaultCohesive()
    self.assertEqual(None, fault.cppHandle)
    return


  def test_adjustTopology(self):
    """
    Test initialize().
    """
    raise NotImplementedError("Need to implement unit test.")
    return
  

# End of file 
