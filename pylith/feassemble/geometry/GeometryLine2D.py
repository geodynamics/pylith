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

## @file pylith/feassemble/geometry/GeometryLine2D.py
##
## @brief Python object for geometry of a 1-D line finite-element cell
## in 2-D.

# ----------------------------------------------------------------------
# GeometryLine2D class
class GeometryLine2D(object):
  """
  Python object for geometry of a 1-D line finite-element cell in 2-D.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.GeometryLine2D()
    return


# End of file 
