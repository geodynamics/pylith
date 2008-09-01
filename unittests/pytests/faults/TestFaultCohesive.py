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


  def test_useFaultMesh(self):
    """
    Test useFaultMesh().
    """
    fault = FaultCohesive()
    fault._configure()
    self.assertEqual(False, fault.useFaultMesh)

    fault.useFaultMesh = True;
    self.assertEqual(True, fault.useFaultMesh)
    return


  def test_faultMeshFilename(self):
    """
    Test faultMeshFilename().
    """
    fault = FaultCohesive()
    fault._configure()
    self.assertEqual("fault.inp", fault.faultMeshFilename)

    filename = "SanAndreas.inp"
    fault.faultMeshFilename = filename
    self.assertEqual(filename, fault.faultMeshFilename)
    return


# End of file 
