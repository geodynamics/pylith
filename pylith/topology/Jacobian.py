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

  def __init__(self, field, matrixType="unknown", blockOkay=False):
    """
    Constructor.

    @param fields Solution fields.
    """
    # If matrix type has not been set, then set it to a value that will work.
    if matrixType == "unknown":
      matrixType = "sbaij"

    #print "MATRIX TYPE: %s, BLOCKOKAY: %s" % (matrixType, blockOkay)
    ModuleJacobian.__init__(self, field, matrixType, blockOkay)
    return
    

  def write(self, filename, comm):
    """
    Write Jacobian to binary file.
    """
    ModuleJacobian.write(self, filename, comm.handle)
    return


  def cleanup(self):
    """
    Dellocate PETSC and local data structures.
    """
    self.deallocate()
    return


# End of file
