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

## @file unittests/pytests/materials/TestElasticStrain1D.py

## @brief Unit testing of ElasticStrain1D object.

import unittest

from pylith.materials.ElasticStrain1D import ElasticStrain1D

# ----------------------------------------------------------------------
class TestElasticStrain1D(unittest.TestCase):
  """
  Unit testing of ElasticStrain1D object.
  """

  def setUp(self):
    """
    Setup test subject.
    """
    self.material = ElasticStrain1D()
    return
  

  def test_constructor(self):
    """
    Test constructor.
    """
    self.assertEqual(1, self.material.dimension())
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
    self.assertEqual(1, self.material.tensorSize())
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
    self.failIf(dt < 1.0e+30)
  

  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.materials.ElasticStrain1D import material
    m = material()
    return


# End of file 
