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

## @file unittests/pytests/topology/TestFieldBase.py

## @brief Unit testing of FieldBase object.

import unittest

from pylith.topology.topology import FieldBase

# ----------------------------------------------------------------------
class TestFieldBase(unittest.TestCase):
  """
  Unit testing of FieldBase object.
  """

  def test_vectorfield(self):
    self.assertEqual(0, FieldBase.SCALAR)
    self.assertEqual(1, FieldBase.VECTOR)
    self.assertEqual(2, FieldBase.TENSOR)
    self.assertEqual(3, FieldBase.OTHER)
    self.assertEqual(4, FieldBase.MULTI_SCALAR)
    self.assertEqual(5, FieldBase.MULTI_VECTOR)
    self.assertEqual(6, FieldBase.MULTI_TENSOR)
    self.assertEqual(7, FieldBase.MULTI_OTHER)
    return


  def test_domain(self):
    self.assertEqual(0, FieldBase.VERTICES_FIELD)
    self.assertEqual(1, FieldBase.CELLS_FIELD)    
    return


# End of file 
