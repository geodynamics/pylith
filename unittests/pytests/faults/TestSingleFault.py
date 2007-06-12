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

## @file unittests/pytests/materials/TestSingleFault.py

## @brief Unit testing of SingleFault object.

import unittest

# ----------------------------------------------------------------------
class TestSingleFault(unittest.TestCase):
  """
  Unit testing of SingleFault object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    from pylith.faults.SingleFault import SingleFault
    faults = SingleFault()
    return


  def test_configure(self):
    """
    Test _configure().
    """
    from pylith.faults.SingleFault import SingleFault
    faults = SingleFault()
    from pylith.faults.Fault import Fault
    faults.inventory.fault = Fault()
    faults._configure()
    self.assertEqual(1, len(faults.bin))
    return


# End of file 
