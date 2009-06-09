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

## @file pylith/topology/Jacobian.py
##
## @brief Python object for system Jacobian.

from topology import Jacobian as ModuleJacobian

# ----------------------------------------------------------------------
# Jacobian class
class Jacobian(ModuleJacobian):
  """
  Python object for system Jacobian.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, fields, matrixType="unknown"):
    """
    Constructor.

    @param fields Solution fields.
    """
    # If matrix type has not been set, then set it to a value that will work.
    if matrixType == "unknown":
      matrixType = "sbaij"

    ModuleJacobian.__init__(self, fields, matrixType)
    return
    

  def cleanup(self):
    """
    Dellocate PETSC and local data structures.
    """
    self.deallocate()
    return


# End of file
