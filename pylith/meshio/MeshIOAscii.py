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

## @file pyre/meshio/MeshIOAscii.py

## @brief Python object for reading/writing finite-element mesh from
## simple ASCII file.

import FIAT.shapes

from pyre.components.Component import Component

# MeshIO class
class MeshIO(Component):
  """
  Python object for reading/writing finite-element mesh from simple
  ASCII file.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="meshio"):
    """Constructor."""
    Component.__init__(self, name, facility="meshio")
    import pylith.meshio.meshio as bindings
    self.cppHandler = bindings.MeshIOAscii()
    return


# version
__id__ = "$Id$"

# End of file 
