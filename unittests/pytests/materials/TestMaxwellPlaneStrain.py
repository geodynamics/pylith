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

## @file unittests/pytests/materials/TestMaxwellPlaneStrain.py

## @brief Unit testing of MaxwellPlaneStrain object.

import unittest

from pylith.materials.MaxwellPlaneStrain import MaxwellPlaneStrain

# ----------------------------------------------------------------------
class TestMaxwellPlaneStrain(unittest.TestCase):
  """
  Unit testing of MaxwellPlaneStrain object.
  """

  def setUp(self):
    """
    Setup test subject.
    """
    self.material = MaxwellPlaneStrain()
    return
  

  def test_constructor(self):
    """
    Test constructor.
    """
    self.assertEqual(2, self.material.dimension())
    return


  def test_useElasticBehavior(self):
    """
    Test useElasticBehavior().
    """
    self.material.useElasticBehavior(False)
    return


  def testHasStateVars(self):
    self.failUnless(self.material.hasStateVars())
    return


  def testTensorSize(self):
    self.assertEqual(3, self.material.tensorSize())
    return


  def testNeedNewJacobian(self):
    """
    Test needNewJacobian().
    """
    # Default should be False.
    self.failIf(self.material.needNewJacobian())

    # Changing time step should require new Jacobian.
    self.material.timeStep(1.0)
    self.material.timeStep(2.0)
    self.failUnless(self.material.needNewJacobian())
    return
  

  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.materials.MaxwellPlaneStrain import material
    m = material()
    return


# End of file 
