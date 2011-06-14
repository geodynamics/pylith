#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2011 University of California, Davis
#
# See COPYING for license information.
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
    self.filter = True
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_vertex_filter():
  """
  Factory associated with VertexFilter.
  """
  return VertexFilterVecNormSubMesh()


# End of file 
