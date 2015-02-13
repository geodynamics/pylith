#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/feassemble/geometry/CellGeometry.py
##
## @brief Python abstract base class for geometry of a finite-element cell.

import feassemble

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
# GeometryLine2D class
class GeometryLine2D(feassemble.GeometryLine2D):
  """
  Python object for geometry of a 1-D finite-element cell in 2-D.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    feassemble.GeometryLine2D.__init__(self)
    return


# ----------------------------------------------------------------------
# GeometryLine3D class
class GeometryLine3D(feassemble.GeometryLine3D):
  """
  Python object for geometry of a 1-D finite-element cell in 3-D.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    feassemble.GeometryLine3D.__init__(self)
    return


# ----------------------------------------------------------------------
# GeometryTri2D class
class GeometryTri2D(feassemble.GeometryTri2D):
  """
  Python object for geometry of a 2-D triangular finite-element cell
  in 2-D.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    feassemble.GeometryTri2D.__init__(self)
    return


# ----------------------------------------------------------------------
# GeometryTri3D class
class GeometryTri3D(feassemble.GeometryTri3D):
  """
  Python object for geometry of a 2-D triangular finite-element cell
  in 3-D.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    feassemble.GeometryTri3D.__init__(self)
    return


# ----------------------------------------------------------------------
# GeometryQuad2D class
class GeometryQuad2D(feassemble.GeometryQuad2D):
  """
  Python object for geometry of a 2-D quadrilateral finite-element cell
  in 2-D.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    feassemble.GeometryQuad2D.__init__(self)
    return


# ----------------------------------------------------------------------
# GeometryQuad3D class
class GeometryQuad3D(feassemble.GeometryQuad3D):
  """
  Python object for geometry of a 2-D quadrilateral finite-element cell
  in 3-D.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    feassemble.GeometryQuad3D.__init__(self)
    return


# ----------------------------------------------------------------------
# GeometryTet3D class
class GeometryTet3D(feassemble.GeometryTet3D):
  """
  Python object for geometry of a 3-D tetrahedral finite-element cell
  in 3-D.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    feassemble.GeometryTet3D.__init__(self)
    return


# ----------------------------------------------------------------------
# GeometryHex3D class
class GeometryHex3D(feassemble.GeometryHex3D):
  """
  Python object for geometry of a 3-D hexahedral finite-element cell
  in 3-D.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    feassemble.GeometryHex3D.__init__(self)
    return


# End of file 
