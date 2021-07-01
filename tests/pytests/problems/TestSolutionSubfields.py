#!/usr/bin/env python
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ======================================================================
#

## @file tests/pytests/problems/TestSolutionSubfields.py

## @brief Unit testing of containers with solution subfields.

import unittest


# ----------------------------------------------------------------------
class TestSolutionSubfields(unittest.TestCase):
  """Unit testing of containers with solution subfields.
  """

  def test_SolnDisp(self):
    """Test SolnDisp.
    """
    from pylith.problems.SolnDisp import SolnDisp
    subfields = SolnDisp()
    subfields._configure()

    components = subfields.components()
    self.assertEqual(subfields.displacement, components[0])
    return
    
  
  def test_SolnDispLagrange(self):
    """Test SolnDispLagrange.
    """
    from pylith.problems.SolnDispLagrange import SolnDispLagrange
    subfields = SolnDispLagrange()
    subfields._configure()

    components = subfields.components()
    self.assertEqual(subfields.displacement, components[0])
    self.assertEqual(subfields.lagrangeFault, components[1])
    return
    
  
  def test_SolnDispVel(self):
    """Test SolnDispVel.
    """
    from pylith.problems.SolnDispVel import SolnDispVel
    subfields = SolnDispVel()
    subfields._configure()

    components = subfields.components()
    self.assertEqual(subfields.displacement, components[0])
    self.assertEqual(subfields.velocity, components[1])
    return
    
  
  def test_SolnDispVelLagrange(self):
    """Test SolnDispVel.
    """
    from pylith.problems.SolnDispVelLagrange import SolnDispVelLagrange
    subfields = SolnDispVelLagrange()
    subfields._configure()

    components = subfields.components()
    self.assertEqual(subfields.displacement, components[0])
    self.assertEqual(subfields.velocity, components[1])
    self.assertEqual(subfields.lagrangeFault, components[2])
    return
    
  
  def test_SolnDispPres(self):
    """Test SolnDispPres.
    """
    from pylith.problems.SolnDispPres import SolnDispPres
    subfields = SolnDispPres()
    subfields._configure()

    components = subfields.components()
    self.assertEqual(subfields.displacement, components[0])
    self.assertEqual(subfields.pressure, components[1])
    return
    
  
  def test_SolnDispPresLagrange(self):
    """Test SolnDispPres.
    """
    from pylith.problems.SolnDispPresLagrange import SolnDispPresLagrange
    subfields = SolnDispPresLagrange()
    subfields._configure()

    components = subfields.components()
    self.assertEqual(subfields.displacement, components[0])
    self.assertEqual(subfields.pressure, components[1])
    self.assertEqual(subfields.lagrangeFault, components[2])
    return
    
  
# End of file 
