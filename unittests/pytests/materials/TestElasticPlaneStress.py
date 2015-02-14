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

## @file unittests/pytests/materials/TestElasticPlaneStress.py

## @brief Unit testing of ElasticPlaneStress object.

import unittest

from pylith.materials.ElasticPlaneStress import ElasticPlaneStress

# ----------------------------------------------------------------------
class TestElasticPlaneStress(unittest.TestCase):
  """
  Unit testing of ElasticPlaneStress object.
  """

  def setUp(self):
    """
    Setup test subject.
    """
    self.material = ElasticPlaneStress()
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
    self.failIf(self.material.hasStateVars())
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

    # Changing time step should not require new Jacobian.
    self.material.timeStep(1.0)
    self.material.timeStep(2.0)
    self.failIf(self.material.needNewJacobian())
    return
  

  def testStableTimeStepImplicit(self):
    """
    Test stableTimeStepImplicit().
    """
    from pylith.topology.Mesh import Mesh
    mesh = Mesh()
    dt = self.material.stableTimeStepImplicit(mesh)
    from pylith.utils.utils import maxdouble
    self.assertAlmostEqual(1.0, dt/maxdouble())
  

  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.materials.ElasticPlaneStress import material
    m = material()
    return


# End of file 
