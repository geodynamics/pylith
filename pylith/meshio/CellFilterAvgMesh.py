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

## @file pyre/meshio/CellFilterAvgMesh.py
##
## @brief Python class for averageing cell fields over each cell's
## quadrature points when writing finite-element data.
##
## Factory: output_cell_filter

from CellFilter import CellFilter
from meshio import MeshCellFilterAvg as ModuleCellFilterAvg

# CellFilterAvgMesh class
class CellFilterAvgMesh(CellFilter, ModuleCellFilterAvg):
  """
  Python class for average cell fields over each cell's quadrature
  points when writing finite-element data.

  Factory: output_cell_filter
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="cellfilteravgmesh"):
    """
    Constructor.
    """
    CellFilter.__init__(self, name)
    ModuleCellFilterAvg.__init__(self)
    return


  def initialize(self, quadrature):
    """
    Initialize output manager.
    """
    self.quadrature(quadrature)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_cell_filter():
  """
  Factory associated with CellFilter.
  """
  return CellFilterAvgMesh()


# End of file 
