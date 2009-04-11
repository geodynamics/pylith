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

## @file pyre/meshio/VertexFilterVecNorm.py
##
## @brief Python class for computing vector norm for each vertex for
## field over vertices when writing finite-element data.
##
## Factory: output_vertex_filter

from VertexFilter import VertexFilter
from meshio import MeshVertexFilterVecNorm as ModuleMeshObject
from meshio import SubMeshVertexFilterVecNorm as ModuleSubMeshObject

# MeshVertexFilterVecNorm class
class MeshVertexFilterVecNorm(VertexFilter, ModuleMeshObject):
  """
  Python class for computing vector norm for each vertex for field
  over vertices when writing finite-element data.

  Factory: output_vertex_filter
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="meshvertexfiltervecnorm"):
    """
    Constructor.
    """
    VertexFilter.__init__(self, name)
    ModuleMeshObject.__init__(self)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def mesh_output_vertex_filter():
  """
  Factory associated with VertexFilter.
  """
  return MeshVertexFilterVecNorm()


# SubMeshVertexFilterVecNorm class
class SubMeshVertexFilterVecNorm(VertexFilter, ModuleSubMeshObject):
  """
  Python class for computing vector norm for each vertex for field
  over vertices when writing finite-element data.

  Factory: output_vertex_filter
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="meshvertexfiltervecnorm"):
    """
    Constructor.
    """
    VertexFilter.__init__(self, name)
    ModuleSubMeshObject.__init__(self)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def mesh_output_vertex_filter():
  """
  Factory associated with VertexFilter.
  """
  return MeshVertexFilterVecNorm()


def submesh_output_vertex_filter():
  """
  Factory associated with VertexFilter.
  """
  return SubMeshVertexFilterVecNorm()


# End of file 
