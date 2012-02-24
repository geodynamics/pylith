#!/usr/bin/env python
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2012 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file unittests/pytests/meshio/TestDataWriterHDF5Ext.py

## @brief Unit testing of Python DataWriterHDF5Ext object.

import unittest

from pylith.meshio.DataWriterHDF5ExtMesh import DataWriterHDF5ExtMesh
from pylith.meshio.DataWriterHDF5ExtSubMesh import DataWriterHDF5ExtSubMesh
from pylith.meshio.DataWriterHDF5ExtSubSubMesh import DataWriterHDF5ExtSubSubMesh

# ----------------------------------------------------------------------
class TestDataWriterHDF5ExtMesh(unittest.TestCase):
  """
  Unit testing of Python DataWriterHDF5ExtMesh object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    filter = DataWriterHDF5ExtMesh()
    filter._configure()
    return


  def test_initialize(self):
    """
    Test constructor.
    """
    filter = DataWriterHDF5ExtMesh()
    filter._configure()

    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    filter.initialize(normalizer)
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.DataWriterHDF5ExtMesh import data_writer
    filter = data_writer()
    return


# ----------------------------------------------------------------------
class TestDataWriterHDF5ExtSubMesh(unittest.TestCase):
  """
  Unit testing of Python DataWriterHDF5Ext object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    filter = DataWriterHDF5ExtSubMesh()
    filter._configure()
    return


  def test_initialize(self):
    """
    Test constructor.
    """
    filter = DataWriterHDF5ExtSubMesh()
    filter._configure()

    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    filter.initialize(normalizer)
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.DataWriterHDF5ExtSubMesh import data_writer
    filter = data_writer()
    return


# ----------------------------------------------------------------------
class TestDataWriterHDF5ExtSubSubMesh(unittest.TestCase):
  """
  Unit testing of Python DataWriterHDF5Ext object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    filter = DataWriterHDF5ExtSubSubMesh()
    filter._configure()
    return


  def test_initialize(self):
    """
    Test constructor.
    """
    filter = DataWriterHDF5ExtSubSubMesh()
    filter._configure()

    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    filter.initialize(normalizer)
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.DataWriterHDF5ExtSubSubMesh import data_writer
    filter = data_writer()
    return


# End of file 
