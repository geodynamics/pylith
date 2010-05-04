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

## @file tests/3d/tet4/TestAxialElasticIsotropic.py
##
## @brief Test suite for testing pylith with axial compresison in
## z-direction for 3-D box.

import unittest
import numpy
import tables

def run_pylith():
  """
  Run pylith.
  """
  if not "done" in dir(run_pylith):
    from pylith.PyLithApp import PyLithApp
    app = PyLithApp("axialelasticisotropic")
    app.run()
    run_pylith.done = True
  return


class TestAxialElasticIsotropic(unittest.TestCase):
  """
  Test suite for testing pylith with axial compression in z-direction
  for 3-D box.
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
