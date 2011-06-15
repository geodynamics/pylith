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
# Copyright (c) 2010-2011 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file unittests/pytests/meshio/TestVertexFilterVecNorm.py

## @brief Unit testing of Python VertexFilterVecNorm object.

import unittest

from pylith.meshio.VertexFilterVecNormMesh import VertexFilterVecNormMesh
from pylith.meshio.VertexFilterVecNormSubMesh import VertexFilterVecNormSubMesh

# ----------------------------------------------------------------------
class TestVertexFilterVecNormMesh(unittest.TestCase):
  """
  Unit testing of Python VertexFilterVecNorm object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    filter = VertexFilterVecNormMesh()
    filter._configure()
    self.failIf(filter.filter is None)
    return


  def test_initialize(self):
    """
    Test constructor.
    """
    filter = VertexFilterVecNormMesh()
    filter._configure()
    filter.initialize()
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.VertexFilterVecNormMesh import output_vertex_filter
    filter = output_vertex_filter()
    return


# ----------------------------------------------------------------------
class TestVertexFilterVecNormSubMesh(unittest.TestCase):
  """
  Unit testing of Python VertexFilterVecNorm object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    filter = VertexFilterVecNormSubMesh()
    filter._configure()
    return


  def test_initialize(self):
    """
    Test constructor.
    """
    filter = VertexFilterVecNormSubMesh()
    filter._configure()
    filter.initialize()
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.VertexFilterVecNormSubMesh import output_vertex_filter
    filter = output_vertex_filter()
    return


# End of file 
