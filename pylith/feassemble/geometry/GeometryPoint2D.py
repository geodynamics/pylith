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

## @file pylith/feassemble/geometry/GeometryPoint2D.py
##
## @brief Python object for geometry of a 0-D finite-element cell in 2-D.

# ----------------------------------------------------------------------
# GeometryPoint2D class
class GeometryPoint2D(object):
  """
  Python object for geometry of a 0-D finite-element cell in 2-D.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.GeometryPoint2D()
    return


# End of file 
