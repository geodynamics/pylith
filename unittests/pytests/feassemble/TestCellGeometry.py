#!/usr/bin/env python
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file unittests/pytests/feassemble/TestCellGeometry.py

## @brief Unit testing of Python CellGeometry object.

import unittest
import numpy

from pylith.feassemble.CellGeometry import GeometryLine2D


# ----------------------------------------------------------------------
class TestCellGeometry(unittest.TestCase):
  """
  Unit testing of Python CellGeometry object.
  """

  def test_constructors(self):
    """
    Test constructors for cell geometry objects.
    """
    
    geometry = GeometryLine2D()
    self.assertEqual(1, geometry.cellDim())
    self.assertEqual(2, geometry.spaceDim())

    return
    

# End of file 
