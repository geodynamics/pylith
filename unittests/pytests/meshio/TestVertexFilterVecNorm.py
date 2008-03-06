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
    return


  def test_initialize(self):
    """
    Test constructor.
    """
    filter = VertexFilterVecNorm()
    filter._configure()
    filter.initialize()
    self.assertNotEqual(None, filter.cppHandle)    
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.VertexFilterVecNorm import output_vertex_filter
    filter = output_vertex_filter()
    return


# End of file 
