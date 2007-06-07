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

## @file pylith/feassemble/geometry/GeometryTet3D.py
##
## @brief Python object for geometry of a 3-D tetrahedral
## finite-element cell in 3-D.

# ----------------------------------------------------------------------
# GeometryTet3D class
class GeometryTet3D(object):
  """
  Python object for geometry of a 3-D tetrahedral finite-element cell
  in 3-D.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.GeometryTet3D()
    return


# End of file 
