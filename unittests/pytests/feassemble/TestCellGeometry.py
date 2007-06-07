#!/usr/bin/env python
#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# {LicenseText}
#
# ======================================================================
#

## @file unittests/pytests/feassemble/TestCellGeometry.py

## @brief Unit testing of Python CellGeometry object.

import unittest
import numpy

from pylith.feassemble.geometry.GeometryPoint1D import GeometryPoint1D
from pylith.feassemble.geometry.GeometryPoint2D import GeometryPoint2D
from pylith.feassemble.geometry.GeometryPoint3D import GeometryPoint3D
from pylith.feassemble.geometry.GeometryLine1D import GeometryLine1D
from pylith.feassemble.geometry.GeometryLine2D import GeometryLine2D
from pylith.feassemble.geometry.GeometryLine3D import GeometryLine3D
from pylith.feassemble.geometry.GeometryTri2D import GeometryTri2D
from pylith.feassemble.geometry.GeometryTri3D import GeometryTri3D
from pylith.feassemble.geometry.GeometryQuad2D import GeometryQuad2D
from pylith.feassemble.geometry.GeometryQuad3D import GeometryQuad3D
from pylith.feassemble.geometry.GeometryTet3D import GeometryTet3D
from pylith.feassemble.geometry.GeometryHex3D import GeometryHex3D


# ----------------------------------------------------------------------
class TestCellGeometry(unittest.TestCase):
  """
  Unit testing of Python CellGeometry object.
  """

  def test_constructors(self):
    """
    Test constructors for cell geometry objects.
    """
    
    geometry = GeometryPoint1D()
    self.failIfEqual(None, geometry.cppHandle)

    geometry = GeometryPoint2D()
    self.failIfEqual(None, geometry.cppHandle)

    geometry = GeometryPoint3D()
    self.failIfEqual(None, geometry.cppHandle)

    geometry = GeometryLine1D()
    self.failIfEqual(None, geometry.cppHandle)

    geometry = GeometryLine2D()
    self.failIfEqual(None, geometry.cppHandle)

    geometry = GeometryLine3D()
    self.failIfEqual(None, geometry.cppHandle)

    geometry = GeometryTri2D()
    self.failIfEqual(None, geometry.cppHandle)

    geometry = GeometryTri3D()
    self.failIfEqual(None, geometry.cppHandle)

    geometry = GeometryQuad2D()
    self.failIfEqual(None, geometry.cppHandle)

    geometry = GeometryQuad3D()
    self.failIfEqual(None, geometry.cppHandle)

    geometry = GeometryTet3D()
    self.failIfEqual(None, geometry.cppHandle)

    geometry = GeometryHex3D()
    self.failIfEqual(None, geometry.cppHandle)

    return
    

# End of file 
