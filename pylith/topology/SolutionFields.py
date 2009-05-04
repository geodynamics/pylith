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

## @file pylith/topology/SolutionFields.py
##
## @brief Python object for managing fields associated with problem
## solution.

from topology import SolutionFields as ModuleSolutionFields

# ----------------------------------------------------------------------
# SolutionFields class
class SolutionFields(ModuleSolutionFields):
  """
  Python object for managing fields associated with problem solution.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, mesh):
    """
    Constructor.
    """
    ModuleSolutionFields.__init__(self, mesh)
    return
    

  def cleanup(self):
    """
    Deallocate PETSc and local data structures.
    """
    self.deallocate()
    return
    

# End of file
