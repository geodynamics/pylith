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

## @file pyre/meshio/VertexFilter.py
##
## @brief Python abstract base class for filtering vertex fields when
## writing finite-element data.
##
## Factory: output_vertex_filter

from pylith.utils.PetscComponent import PetscComponent

# VertexFilter class
class VertexFilter(PetscComponent):
  """
  Python abstract base class for filtering cell fields when writing
  finite-element data.

  Factory: output_vertex_filter
  """

  # INVENTORY //////////////////////////////////////////////////////////

  # None

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="vertexfilter"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="vertexfilter")
    self.filter = None
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    return


  def initialize(self):
    """
    Initialize output manager.
    """
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_vertex_filter():
  """
  Factory associated with VertexFilter.
  """
  return VertexFilter()


# End of file 
