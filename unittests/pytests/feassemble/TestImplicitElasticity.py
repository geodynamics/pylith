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

## @file unittests/pytests/feassemble/TestImplicitElasticity.py

## @brief Unit testing of Python ImplicitElasticity object.

import unittest
from pylith.feassemble.ImplicitElasticity import ImplicitElasticity

# ----------------------------------------------------------------------
class TestImplicitElasticity(unittest.TestCase):
  """
  Unit testing of Python ImplicitElasticity object.
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

    integrator = ImplicitElasticity()
    integrator.initQuadrature(q)
    self.assertEqual(minJacobian, integrator.quadrature.minJacobian)
    return
    

# End of file 
