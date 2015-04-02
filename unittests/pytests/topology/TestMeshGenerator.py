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

## @file unittests/pytests/topology/TestMeshGenerator.py
##
## @brief Unit testing of MeshGenerator object.

import unittest

from pylith.topology.MeshGenerator import MeshGenerator

# ----------------------------------------------------------------------
class TestMeshGenerator(unittest.TestCase):
  """
  Unit testing of MeshGenerator object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    generator = MeshGenerator()
    return
  

  def test_debug(self):
    """
    Test debug().
    """
    generator = MeshGenerator()

    value = False # default should be False
    self.failUnless(value == generator.debug)

    value = True
    generator.debug = value
    self.failUnless(value == generator.debug)
    return


  def test_interpolate(self):
    """
    Test interpolate access.
    """
    generator = MeshGenerator()

    value = False # default should be False
    self.assertEqual(value, generator.interpolate)

    value = True
    generator.interpolate = value
    self.assertEqual(value, generator.interpolate)
    return


# End of file 
