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

## @file pylith/feassemble/geometry/GeometryLine3D.py
##
## @brief Python object for geometry of a 1-D line finite-element cell in 3-D.

# ----------------------------------------------------------------------
# GeometryLine3D class
class GeometryLine3D(object):
  """
  Python object for geometry of a 1-D finite-element cell in 3-D.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.GeometryLine3D()
    return


# End of file 
