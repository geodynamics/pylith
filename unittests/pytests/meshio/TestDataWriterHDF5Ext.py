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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file unittests/pytests/meshio/TestDataWriterHDF5Ext.py

## @brief Unit testing of Python DataWriterHDF5Ext object.

import unittest

from pylith.meshio.DataWriterHDF5Ext import DataWriterHDF5Ext

# ----------------------------------------------------------------------
class TestDataWriterHDF5Ext(unittest.TestCase):
  """
  Unit testing of Python DataWriterHDF5Ext object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    filter = DataWriterHDF5Ext()
    filter._configure()
    return


  def test_initialize(self):
    """
    Test constructor.
    """
    filter = DataWriterHDF5Ext()
    filter._configure()

    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    filter.initialize(normalizer)
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.DataWriterHDF5Ext import data_writer
    filter = data_writer()
    return


# End of file 
