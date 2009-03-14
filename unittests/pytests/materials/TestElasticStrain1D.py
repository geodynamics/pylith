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


  def testStableTimeStepImplicit(self):
    maxfloat = 1.0e+30
    self.assertEqual(maxfloat, self.material.stableTimeStepImplicit())
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
  

# End of file 
