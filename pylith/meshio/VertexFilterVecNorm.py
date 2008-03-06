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

# VertexFilterVecNorm class
class VertexFilterVecNorm(VertexFilter):
  """
  Python class for computing vector norm for each vertex for field
  over vertices when writing finite-element data.

  Factory: output_vertex_filter
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="vertexfiltervecnorm"):
    """
    Constructor.
    """
    VertexFilter.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _createCppHandle(self):
    """
    Create handle to C++ object.
    """
    if None == self.cppHandle:
      import meshio as bindings
      self.cppHandle = bindings.VertexFilterVecNorm()
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def output_vertex_filter():
  """
  Factory associated with VertexFilter.
  """
  return VertexFilterVecNorm()


# End of file 
