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
    return


  def test_initialize(self):
    """
    Test constructor.
    """
    from pylith.feassemble.quadrature.Quadrature1D import Quadrature1D
    
    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATSimplex()
    cell.shape = "line"
    cell.degree = 2
    cell.order = 2

    quadrature = Quadrature1D()
    quadrature.cell = cell
    quadrature.preinitialize()
    quadrature.initialize()

    filter = CellFilterAvg()
    filter._configure()
    filter.initialize(quadrature)
    self.assertNotEqual(None, filter.cppHandle)    
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.CellFilterAvg import output_cell_filter
    filter = output_cell_filter()
    return


# End of file 
