#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

## @file tests/2d/quad4/TestDislocation.py
##
## @brief Test suite for testing pylith with shear dislocation for 2-D
## box.

import unittest
import numpy
import tables

def run_pylith():
  """
  Run pylith.
  """
  if not "done" in dir(run_pylith):
    from pylith.PyLithApp import PyLithApp
    app = PyLithApp("shearplanestrain")
    app.run()
    run_pylith.done = True
  return


class TestDislocation(unittest.TestCase):
  """
  Test suite for testing pylith with shear dislocation for 2-D box.
  """

  def setUp(self):
    """
    Setup for test.
    """
    run_pylith()
    return


  def test_disp(self):
    """
    Check displacement field.
    """
    return


# End of file 
