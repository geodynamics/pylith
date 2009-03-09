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

## @file pylith/feassemble/geometry/CellGeometry.py
##
## @brief Python abstract base class for geometry of a finite-element cell.

# ----------------------------------------------------------------------
# CellGeometry class
class CellGeometry(object):
  """
  Python abstract base class for geometry of a finite-element cell.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    return


# ----------------------------------------------------------------------
from feassemble import GeometryPoint1D as ModuleGeometryPoint1D

# GeometryPoint1D class
class GeometryPoint1D(ModuleGeometryPoint1D):
  """
  Python object for geometry of a 0-D finite-element cell in 1-D.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    ModuleGeometryPoint1D.__init__(self)
    return


# End of file 
