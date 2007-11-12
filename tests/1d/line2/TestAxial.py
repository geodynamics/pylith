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

## @file tests/1d/line2/TestAxial.py
##
## @brief Test suite for testing pylith with 1-D axial extension.

import unittest
#import numpy
#import tables


from pylith.PyLithApp import PyLithApp
class AxialExtension(PyLithApp):

  def __init__(self):
    PyLithApp.__init__(self, "axialextension")
    return


def run_pylith():
  """
  Run pylith.
  """
  if not "done" in dir(run_pylith):
    app = AxialExtension()
    app.run()
    run_pylith.done = True
  return


class TestAxial(unittest.TestCase):
  """
  Test suite for testing pylith with 1-D axial extension.
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
