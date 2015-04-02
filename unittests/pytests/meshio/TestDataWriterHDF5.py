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

## @file unittests/pytests/meshio/TestDataWriterHDF5.py

## @brief Unit testing of Python DataWriterHDF5 object.

import unittest

from pylith.meshio.DataWriterHDF5 import DataWriterHDF5

# ----------------------------------------------------------------------
class TestDataWriterHDF5(unittest.TestCase):
  """
  Unit testing of Python DataWriterHDF5 object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    filter = DataWriterHDF5()
    filter._configure()
    return


  def test_initialize(self):
    """
    Test constructor.
    """
    filter = DataWriterHDF5()
    filter._configure()

    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    filter.initialize(normalizer)
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.DataWriterHDF5 import data_writer
    filter = data_writer()
    return


# End of file 
