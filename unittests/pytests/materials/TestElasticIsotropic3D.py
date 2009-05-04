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

## @file unittests/pytests/materials/TestElasticIsotropic3D.py

## @brief Unit testing of ElasticIsotropic3D object.

import unittest

from pylith.materials.ElasticIsotropic3D import ElasticIsotropic3D

# ----------------------------------------------------------------------
class TestElasticIsotropic3D(unittest.TestCase):
  """
  Unit testing of ElasticIsotropic3D object.
  """

  def setUp(self):
    """
    Setup test subject.
    """
    self.material = ElasticIsotropic3D()
    return
  

  def test_constructor(self):
    """
    Test constructor.
    """
    self.assertEqual(3, self.material.dimension())
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
    self.assertEqual(6, self.material.tensorSize())
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
