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

## @file pylith/feassemble/geometry/GeometryPoint3D.py
##
## @brief Python object for geometry of a 0-D finite-element cell in 3-D.

# ----------------------------------------------------------------------
# GeometryPoint3D class
class GeometryPoint3D(object):
  """
  Python object for geometry of a 0-D finite-element cell in 3-D.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.GeometryPoint3D()
    return


# End of file 
