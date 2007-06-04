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

## @file unittests/pytests/feassemble/TestIntegrator.py

## @brief Unit testing of Python Integrator object.

import unittest
from pylith.feassemble.Integrator import Integrator

# ----------------------------------------------------------------------
class TestIntegrator(unittest.TestCase):
  """
  Unit testing of Python Integrator object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    i = Integrator()
    return


  def test_finalize(self):
    """
    Test constructor.
    """
    i = Integrator()
    i.finalize()
    return


# End of file 
