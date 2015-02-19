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

## @file unittests/pytests/materials/TestHomogeneous.py

## @brief Unit testing of Homogenous object.

import unittest

# ----------------------------------------------------------------------
class TestHomogeneous(unittest.TestCase):
  """
  Unit testing of Homogeneous object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    from pylith.materials.Homogeneous import Homogeneous
    materials = Homogeneous()
    return


  def test_configure(self):
    """
    Test _configure().
    """
    from pylith.materials.Homogeneous import Homogeneous
    materials = Homogeneous()
    from pylith.materials.ElasticIsotropic3D import ElasticIsotropic3D
    materials.inventory.material = ElasticIsotropic3D()
    materials._configure()
    self.assertEqual(1, len(materials.components()))
    return


# End of file 
