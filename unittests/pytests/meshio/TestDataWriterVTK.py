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

## @file unittests/pytests/meshio/TestDataWriterVTK.py

## @brief Unit testing of Python DataWriterVTK object.

import unittest

from pylith.meshio.DataWriterVTK import DataWriterVTK

# ----------------------------------------------------------------------
class TestDataWriterVTK(unittest.TestCase):
  """
  Unit testing of Python DataWriterVTK object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    filter = DataWriterVTK()
    filter._configure()
    return


  def test_initialize(self):
    """
    Test constructor.
    """
    filter = DataWriterVTK()
    filter._configure()

    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    filter.initialize(normalizer)
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.DataWriterVTK import data_writer
    filter = data_writer()
    return


# End of file 
