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

## @file pylith/topology/ReverseCuthillMcKee.py
##
## @brief Python interface to PETSc reverse Cuthill-McKee reordering
## of mesh cells and vertices.

from topology import ReverseCuthillMcKee as ModuleReverseCuthillMcKee

# ReverseCuthillMcKee class
class ReverseCuthillMcKee(ModuleReverseCuthillMcKee):
  """
  Python interface to PETSc reverse Cuthill-McKee reordering of mesh
  cells and vertices.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    return


  def reorder(self, mesh):
    """
    Set communicator.
    """
    ModuleReverseCuthillMcKee.reorder(mesh)
    return


# End of file
