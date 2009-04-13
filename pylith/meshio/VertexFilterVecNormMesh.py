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

## @file pyre/meshio/VertexFilterVecNormMesh.py
##
## @brief Python class for computing vector norm for each vertex for
## field over vertices when writing finite-element data.
##
## Factory: output_vertex_filter

from VertexFilter import VertexFilter
from meshio import MeshVertexFilterVecNorm as ModuleVertexFilterVecNorm

# VertexFilterVecNormMesh class
class VertexFilterVecNormMesh(VertexFilter, ModuleVertexFilterVecNorm):
  """
  Python class for computing vector norm for each vertex for field
  over vertices when writing finite-element data.

  Factory: output_vertex_filter
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="vertexfiltervecnormmesh"):
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
  return VertexFilterVecNormMesh()


# End of file 
