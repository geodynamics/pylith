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

## @file tests/2d/tri3/TestAxialPlaneStrain.py
##
## @brief Test suite for testing pylith with axial extension in
## y-direction for 2-D box.

import unittest

# Local version of PyLithApp
from pylith.PyLithApp import PyLithApp
class AxialPlaneStrainApp(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="axialplanestrain")
    return


# Helper function to run PyLith
def run_pylith():
  """
  Run pylith.
  """
  if not "done" in dir(run_pylith):
    app = AxialPlaneStrainApp()
    app.run()
    run_pylith.done = True
  return


class TestAxialPlaneStrain(unittest.TestCase):
  """
  Test suite for testing pylith with axial extension in y-direction
  for 2-D box.
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
