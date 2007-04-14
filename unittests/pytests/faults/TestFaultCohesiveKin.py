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

from pylith.faults.FaultCohesiveKin import FaultCohesiveKin

# ----------------------------------------------------------------------
class TestFaultCohesiveKin(unittest.TestCase):
  """
  Unit testing of Fault object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    fault = FaultCohesiveKin()
    self.failIfEqual(None, fault.cppHandle)
    return


# End of file 
