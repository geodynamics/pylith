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

## @file unittests/pytests/meshio/TestCellFilterAvg.py

## @brief Unit testing of Python CellFilterAvg object.

import unittest

from pylith.meshio.CellFilterAvg import CellFilterAvg

# ----------------------------------------------------------------------
class TestCellFilterAvg(unittest.TestCase):
  """
  Unit testing of Python CellFilterAvg object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    filter = CellFilterAvg()
    filter._configure()
    self.failIf(filter.filter is None)
    return


  def test_initialize(self):
    """
    Test constructor.
    """
    spaceDim = 2
    
    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATSimplex()
    cell.inventory.dimension = 2
    cell.inventory.degree = 2
    cell.inventory.order = 2
    cell._configure()

    from pylith.feassemble.Quadrature import Quadrature
    quadrature = Quadrature()
    quadrature.inventory.cell = cell
    quadrature._configure()
    quadrature.preinitialize(spaceDim)
    quadrature.initialize()

    filter = CellFilterAvg()
    filter._configure()
    filter.initialize(quadrature)
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.CellFilterAvg import output_cell_filter
    filter = output_cell_filter()
    return


# End of file 
