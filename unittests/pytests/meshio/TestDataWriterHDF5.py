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
# Copyright (c) 2010-2013 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file unittests/pytests/meshio/TestDataWriterHDF5.py

## @brief Unit testing of Python DataWriterHDF5 object.

import unittest

from pylith.meshio.DataWriterHDF5Mesh import DataWriterHDF5Mesh

# ----------------------------------------------------------------------
class TestDataWriterHDF5Mesh(unittest.TestCase):
  """
  Unit testing of Python DataWriterHDF5Mesh object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    filter = DataWriterHDF5Mesh()
    filter._configure()
    return


  def test_initialize(self):
    """
    Test constructor.
    """
    filter = DataWriterHDF5Mesh()
    filter._configure()

    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    filter.initialize(normalizer)
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.DataWriterHDF5Mesh import data_writer
    filter = data_writer()
    return


# End of file 
