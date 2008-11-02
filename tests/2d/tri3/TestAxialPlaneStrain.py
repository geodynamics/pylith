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
from pylith.utils.VTKDataReader import has_vtk
from pylith.utils.VTKDataReader import VTKDataReader

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
    if has_vtk():
      self.reader = VTKDataReader()
    else:
      self.reader = None
    return

  
  def test_elastic_info(self):
    """
    Check elastic info.
    """
    if self.reader is None:
      return

    data = self.reader.read("axialplanestrain-statevars-elastic_info.vtk")
    return


  def test_soln(self):
    """
    Check solution (displacement) field.
    """
    return


  def test_elastic_statevars(self):
    """
    Check elastic state variables.
    """
    return


# End of file 
