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

from pyre.components.Component import Component

# EmptyBin class
class EmptyBin(Component):
  """
  Python container for a collection of objects.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="emptybin"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="empty_bin")
    return


# End of file 
