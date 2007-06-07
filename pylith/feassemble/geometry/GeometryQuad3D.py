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

## @file pylith/feassemble/geometry/GeometryQuad3D.py
##
## @brief Python object for geometry of a 2-D quadrilateral
## finite-element cell in 3-D.

# ----------------------------------------------------------------------
# GeometryQuad3D class
class GeometryQuad3D(object):
  """
  Python object for geometry of a 2-D quadrilateral finite-element
  cell in 3-D.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.GeometryQuad3D()
    return


# End of file 
