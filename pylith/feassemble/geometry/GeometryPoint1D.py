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

## @file pylith/feassemble/geometry/GeometryPoint1D.py
##
## @brief Python object for geometry of a 0-D finite-element cell in 1-D.

# ----------------------------------------------------------------------
# GeometryPoint1D class
class GeometryPoint1D(object):
  """
  Python object for geometry of a 0-D finite-element cell in 1-D.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.GeometryPoint1D()
    return


# End of file 
