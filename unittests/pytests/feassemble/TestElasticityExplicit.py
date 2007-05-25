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

## @file unittests/pytests/feassemble/TestElasticityExplicit.py

## @brief Unit testing of Python ElasticityExplicit object.

import unittest
from pylith.feassemble.ElasticityExplicit import ElasticityExplicit

# ----------------------------------------------------------------------
class TestElasticityExplicit(unittest.TestCase):
  """
  Unit testing of Python ElasticityExplicit object.
  """

  def test_initQuadrature(self):
    """
    Test initQuadrature().
    """
    from pylith.feassemble.quadrature.Quadrature2D import Quadrature2D
    q = Quadrature2D()
    minJacobian = 4.0e-02;
    q.minJacobian = minJacobian
    from pylith.feassemble.FIATSimplex import FIATSimplex
    q.cell = FIATSimplex()
    q.cell.shape = "triangle"
    q.cell.order = 1
    q.cell.degree = 1

    integrator = ElasticityExplicit()
    integrator.initQuadrature(q)
    self.assertEqual(minJacobian, integrator.quadrature.minJacobian)
    return
    

  def test_needNewJacobian(self):
    """
    Test needNewJacobian().
    """
    integrator = ElasticityExplicit()
    self.assertEqual(False, integrator.needNewJacobian())
    return


# End of file 
