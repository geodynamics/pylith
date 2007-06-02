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

## @file pylith/utils/ObjectBin.py
##
## @brief Python container for a collection of objects.
##
## Factory: object_bin

from pyre.components.Component import Component

# ObjectBin class
class ObjectBin(Component):
  """
  Python container for a collection of objects.

  Factory: object_bin
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="objectbin"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="object_bin")
    self.bin = []
    return


# End of file 
