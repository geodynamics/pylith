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

## @file pylith/feassemble/geometry/CellGeometry.py
##
## @brief Python abstract base class for geometry of a finite-element cell.

# ----------------------------------------------------------------------
# CellGeometry class
class CellGeometry(object):
  """
  Python abstract base class for geometry of a finite-element cell.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    self.cppHandle = None
    return


# End of file 
