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

## @file unittests/pytests/feassemble/TestExplicitElasticity.py

## @brief Unit testing of Python ExplicitElasticity object.

import unittest
from pylith.feassemble.ExplicitElasticity import ExplicitElasticity

# ----------------------------------------------------------------------
class TestExplicitElasticity(unittest.TestCase):
  """
  Unit testing of Python ExplicitElasticity object.
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

    integrator = ExplicitElasticity()
    integrator.initQuadrature(q)
    self.assertEqual(minJacobian, integrator.quadrature.minJacobian)
    return
    

# End of file 
