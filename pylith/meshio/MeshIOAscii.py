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

from MeshIO import MeshIO

# MeshIOAscii class
class MeshIOAscii(MeshIO):
  """
  Python object for reading/writing finite-element mesh from simple
  ASCII file.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="meshioascii"):
    """Constructor."""
    MeshIO.__init__(self, name)
    import pylith.meshio.meshio as bindings
    self.cppHandler = bindings.MeshIOAscii()
    return


# version
__id__ = "$Id$"

# End of file 
