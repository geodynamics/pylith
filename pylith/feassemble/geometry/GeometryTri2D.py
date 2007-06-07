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

## @file pylith/feassemble/geometry/GeometryTri2D.py
##
## @brief Python object for geometry of a 2-D triangular
## finite-element cell in 2-D.

# ----------------------------------------------------------------------
# GeometryTri2D class
class GeometryTri2D(object):
  """
  Python object for geometry of a 2-D triangular finite-element cell
  in 2-D.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.GeometryTri2D()
    return


# End of file 
