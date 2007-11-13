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

## @file tests/1d/line2/TestDislocation.py
##
## @brief Test suite for testing pylith with static dislocation in 1-D
## mesh.

import unittest
#import numpy
#import tables


from pylith.PyLithApp import PyLithApp
class Dislocation(PyLithApp):

  def __init__(self):
    PyLithApp.__init__(self, "dislocation_static")
    return


def run_pylith():
  """
  Run pylith.
  """
  if not "done" in dir(run_pylith):
    app = Dislocation()
    app.run()
    run_pylith.done = True
  return


class TestDislocationStatic(unittest.TestCase):
  """
  Test suite for testing pylith with static dislocation in 1-D mesh.
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
