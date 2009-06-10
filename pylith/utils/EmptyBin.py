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

## @file pylith/utils/EmptyBin.py
##
## @brief Python container for a collection of objects.
##
## Factory: object_bin

from pylith.utils.PetscComponent import PetscComponent

# EmptyBin class
class EmptyBin(PetscComponent):
  """
  Python container for a collection of objects.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="emptybin"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="empty_bin")
    return


# End of file 
