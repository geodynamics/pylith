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

from pylith.meshio.CellFilterAvg import MeshCellFilterAvg
from pylith.meshio.CellFilterAvg import SubMeshCellFilterAvg

# ----------------------------------------------------------------------
class TestMeshCellFilterAvg(unittest.TestCase):
  """
  Unit testing of Python CellFilterAvg object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    filter = MeshCellFilterAvg()
    filter._configure()
    return


  def test_initialize(self):
    """
    Test constructor.
    """
    spaceDim = 1
    
    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATSimplex()
    cell.inventory.shape = "line"
    cell.inventory.degree = 2
    cell.inventory.order = 2
    cell._configure()

    from pylith.feassemble.Quadrature import MeshQuadrature
    quadrature = MeshQuadrature()
    quadrature.inventory.cell = cell
    quadrature._configure()
    quadrature.preinitialize(spaceDim)
    quadrature.initialize()

    filter = MeshCellFilterAvg()
    filter._configure()
    filter.initialize(quadrature)
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.CellFilterAvg import mesh_output_cell_filter
    filter = mesh_output_cell_filter()
    return


# ----------------------------------------------------------------------
class TestSubMeshCellFilterAvg(unittest.TestCase):
  """
  Unit testing of Python CellFilterAvg object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    filter = SubMeshCellFilterAvg()
    filter._configure()
    return


  def test_initialize(self):
    """
    Test constructor.
    """
    spaceDim = 2
    
    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = FIATSimplex()
    cell.inventory.shape = "line"
    cell.inventory.degree = 2
    cell.inventory.order = 2
    cell._configure()

    from pylith.feassemble.Quadrature import SubMeshQuadrature
    quadrature = SubMeshQuadrature()
    quadrature.inventory.cell = cell
    quadrature._configure()
    quadrature.preinitialize(spaceDim)
    quadrature.initialize()

    filter = SubMeshCellFilterAvg()
    filter._configure()
    filter.initialize(quadrature)
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.CellFilterAvg import mesh_output_cell_filter
    filter = mesh_output_cell_filter()
    return


# End of file 
