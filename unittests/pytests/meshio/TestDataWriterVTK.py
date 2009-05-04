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

from pylith.meshio.DataWriterVTKMesh import DataWriterVTKMesh
from pylith.meshio.DataWriterVTKSubMesh import DataWriterVTKSubMesh
from pylith.meshio.DataWriterVTKSubSubMesh import DataWriterVTKSubSubMesh

# ----------------------------------------------------------------------
class TestDataWriterVTKMesh(unittest.TestCase):
  """
  Unit testing of Python DataWriterVTKMesh object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    filter = DataWriterVTKMesh()
    filter._configure()
    return


  def test_initialize(self):
    """
    Test constructor.
    """
    filter = DataWriterVTKMesh()
    filter._configure()

    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    filter.initialize(normalizer)
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.DataWriterVTKMesh import output_data_writer
    filter = output_data_writer()
    return


# ----------------------------------------------------------------------
class TestDataWriterVTKSubMesh(unittest.TestCase):
  """
  Unit testing of Python DataWriterVTK object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    filter = DataWriterVTKSubMesh()
    filter._configure()
    return


  def test_initialize(self):
    """
    Test constructor.
    """
    filter = DataWriterVTKSubMesh()
    filter._configure()

    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    filter.initialize(normalizer)
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.DataWriterVTKSubMesh import output_data_writer
    filter = output_data_writer()
    return


# ----------------------------------------------------------------------
class TestDataWriterVTKSubSubMesh(unittest.TestCase):
  """
  Unit testing of Python DataWriterVTK object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    filter = DataWriterVTKSubSubMesh()
    filter._configure()
    return


  def test_initialize(self):
    """
    Test constructor.
    """
    filter = DataWriterVTKSubSubMesh()
    filter._configure()

    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    filter.initialize(normalizer)
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.DataWriterVTKSubSubMesh import output_data_writer
    filter = output_data_writer()
    return


# End of file 
