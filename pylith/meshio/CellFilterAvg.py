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

## @file pyre/meshio/CellFilterAvg.py
##
## @brief Python class for averageing cell fields over each cell's
## quadrature points when writing finite-element data.
##
## Factory: output_cell_filter

from CellFilter import CellFilter
from meshio import MeshCellFilterAvg as ModuleMeshObject
from meshio import SubMeshCellFilterAvg as ModuleSubMeshObject

# MeshCellFilterAvg class
class MeshCellFilterAvg(CellFilter, ModuleMeshObject):
  """
  Python class for average cell fields over each cell's quadrature
  points when writing finite-element data.

  Factory: mesh_output_cell_filter
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="meshcellfilteravg"):
    """
    Constructor.
    """
    CellFilter.__init__(self, name)
    ModuleMeshObject.__init__(self)
    return


  def initialize(self, quadrature):
    """
    Initialize output manager.
    """
    self.quadrature(quadrature)
    return


# SubMeshCellFilterAvg class
class SubMeshCellFilterAvg(CellFilter, ModuleSubMeshObject):
  """
  Python class for average cell fields over each cell's quadrature
  points when writing finite-element data.

  Factory: submesh_output_cell_filter
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="submeshcellfilteravg"):
    """
    Constructor.
    """
    CellFilter.__init__(self, name)
    ModuleSubMeshObject.__init__(self)
    return


  def initialize(self, quadrature):
    """
    Initialize output manager.
    """
    self.quadrature(quadrature)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def mesh_output_cell_filter():
  """
  Factory associated with CellFilter.
  """
  return MeshCellFilterAvg()


def submesh_output_cell_filter():
  """
  Factory associated with CellFilter.
  """
  return SubMeshCellFilterAvg()


# End of file 
