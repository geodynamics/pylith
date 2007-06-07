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

## @file pylith/feassemble/geometry/GeometryTri3D.py
##
## @brief Python object for geometry of a 2-D triangular
## finite-element cell in 3-D.

# ----------------------------------------------------------------------
# GeometryTri3D class
class GeometryTri3D(object):
  """
  Python object for geometry of a 2-D triangular finite-element cell
  in 3-D.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.GeometryTri3D()
    return


# End of file 
