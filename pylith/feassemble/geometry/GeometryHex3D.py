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

## @file pylith/feassemble/geometry/GeometryHex3D.py
##
## @brief Python object for geometry of a 3-D hexahedral
## finite-element cell in 3-D.

# ----------------------------------------------------------------------
# GeometryHex3D class
class GeometryHex3D(object):
  """
  Python object for geometry of a 3-D hexahedral finite-element cell
  in 3-D.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.GeometryHex3D()
    return


# End of file 
