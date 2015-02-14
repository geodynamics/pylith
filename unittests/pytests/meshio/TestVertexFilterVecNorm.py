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

## @file unittests/pytests/meshio/TestVertexFilterVecNorm.py

## @brief Unit testing of Python VertexFilterVecNorm object.

import unittest

from pylith.meshio.VertexFilterVecNorm import VertexFilterVecNorm

# ----------------------------------------------------------------------
class TestVertexFilterVecNorm(unittest.TestCase):
  """
  Unit testing of Python VertexFilterVecNorm object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    filter = VertexFilterVecNorm()
    filter._configure()
    self.failIf(filter.filter is None)
    return


  def test_initialize(self):
    """
    Test constructor.
    """
    filter = VertexFilterVecNorm()
    filter._configure()
    filter.initialize()
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.VertexFilterVecNorm import output_vertex_filter
    filter = output_vertex_filter()
    return


# End of file 
