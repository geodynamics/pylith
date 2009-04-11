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

## @file unittests/pytests/meshio/TestDataWriterVTK.py

## @brief Unit testing of Python DataWriterVTK object.

import unittest

from pylith.meshio.DataWriterVTK import MeshDataWriterVTK
from pylith.meshio.DataWriterVTK import SubMeshDataWriterVTK

# ----------------------------------------------------------------------
class TestMeshDataWriterVTK(unittest.TestCase):
  """
  Unit testing of Python DataWriterVTK object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    filter = MeshDataWriterVTK()
    filter._configure()
    return


  def test_initialize(self):
    """
    Test constructor.
    """
    filter = MeshDataWriterVTK()
    filter._configure()

    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    filter.initialize(normalizer)
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.DataWriterVTK import mesh_output_data_writer
    filter = mesh_output_data_writer()
    return


# ----------------------------------------------------------------------
class TestSubMeshDataWriterVTK(unittest.TestCase):
  """
  Unit testing of Python DataWriterVTK object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    filter = SubMeshDataWriterVTK()
    filter._configure()
    return


  def test_initialize(self):
    """
    Test constructor.
    """
    filter = SubMeshDataWriterVTK()
    filter._configure()

    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    filter.initialize(normalizer)
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.DataWriterVTK import submesh_output_data_writer
    filter = submesh_output_data_writer()
    return


# End of file 
