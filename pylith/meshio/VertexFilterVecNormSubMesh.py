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

## @file pyre/meshio/VertexFilterVecNormSubMesh.py
##
## @brief Python class for computing vector norm for each vertex for
## field over vertices when writing finite-element data.
##
## Factory: output_vertex_filter

from VertexFilter import VertexFilter
from meshio import SubMeshVertexFilterVecNorm as ModuleVertexFilterVecNorm

# VertexFilterVecNormSubMesh class
class VertexFilterVecNormSubMesh(VertexFilter, ModuleVertexFilterVecNorm):
  """
  Python class for computing vector norm for each vertex for field
  over vertices when writing finite-element data.

  Factory: output_vertex_filter
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="vertexfiltervecnormsubmesh"):
    """
    Constructor.
    """
    VertexFilter.__init__(self, name)
    ModuleVertexFilterVecNorm.__init__(self)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_vertex_filter():
  """
  Factory associated with VertexFilter.
  """
  return VertexFilterVecNormSubMesh()


# End of file 
