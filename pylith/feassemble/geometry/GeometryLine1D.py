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

## @file pylith/feassemble/geometry/GeometryLine1D.py
##
## @brief Python object for geometry of a 1-D line finite-element cell
## in 1-D.

# ----------------------------------------------------------------------
# GeometryLine1D class
class GeometryLine1D(object):
  """
  Python object for geometry of a 1-D line finite-element cell in 1-D.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.GeometryLine1D()
    return


# End of file 
