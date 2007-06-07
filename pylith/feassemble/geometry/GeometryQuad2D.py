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

## @file pylith/feassemble/geometry/GeometryQuad2D.py
##
## @brief Python object for geometry of a 2-D quadrilateral
## finite-element cell in 2-D.

# ----------------------------------------------------------------------
# GeometryQuad2D class
class GeometryQuad2D(object):
  """
  Python object for geometry of a 2-D quadrilateral finite-element
  cell in 2-D.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.GeometryQuad2D()
    return


# End of file 
