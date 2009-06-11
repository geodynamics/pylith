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

## @file unittests/pytests/materials/TestElasticPlaneStrain.py

## @brief Unit testing of ElasticPlaneStrain object.

import unittest

from pylith.materials.ElasticPlaneStrain import ElasticPlaneStrain

# ----------------------------------------------------------------------
class TestElasticPlaneStrain(unittest.TestCase):
  """
  Unit testing of ElasticPlaneStrain object.
  """

  def setUp(self):
    """
    Setup test subject.
    """
    self.material = ElasticPlaneStrain()
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
    self.failIf(dt < 1.0e+30)
  

  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.materials.ElasticPlaneStrain import material
    m = material()
    return


# End of file 
