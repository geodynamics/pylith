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

## @file pylith/materials/MaterialsBin.py

## @brief Python container for materials.

from pyre.components.Component import Component

# MaterialsBin class
class MaterialsBin(Component):
  """
  Python container for materials.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="materialsbin"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="materialsbin")
    self.materials = []
    return


# End of file 
