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
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file tests/pytests/problems/TestSolution.py

## @brief Unit testing of Solution object.

import unittest


from pylith.problems.Solution import Solution

from pythia.pyre.units.time import second
from pythia.pyre.units.length import km
from pythia.pyre.units.pressure import MPa
timeScale = 2.0*second
lengthScale = 3.0*km
pressureScale = 0.2*MPa
spaceDim = 2

# ----------------------------------------------------------------------
class TestSolution(unittest.TestCase):
  """
  Unit testing of Solution object.
  """

  def setUp(self):
    from pylith.problems.SolnDispVelLagrange import SolnDispVelLagrange
    solution = Solution()
    solution.inventory.subfields = SolnDispVelLagrange()
    solution.inventory.subfields.inventory.displacement._configure()
    solution.inventory.subfields.inventory.velocity._configure()
    solution.inventory.subfields.inventory.lagrangeFault._configure()
    solution.inventory.subfields._configure()
    solution._configure()
    self.solution = solution
    return
  

  def test_constructor(self):
    """
    Test constructor.
    """
    solution = Solution()
    solution._configure()
    return
    
  
  def test_preinitialize(self):
    """
    Test preinitialize().
    """
    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    normalizer._configure()
    normalizer.setTimeScale(timeScale)
    normalizer.setLengthScale(lengthScale)
    normalizer.setPressureScale(pressureScale)

    self.solution.preinitialize(normalizer, spaceDim)

    # Displacement
    from pylith.topology.Field import Field
    self.assertEqual(spaceDim, self.solution.subfields.displacement.ncomponents)
    self.assertEqual(Field.VECTOR, self.solution.subfields.displacement.vectorFieldType)
    self.assertEqual(lengthScale, self.solution.subfields.displacement.scale)

    # Velocity
    from pylith.topology.Field import Field
    self.assertEqual(spaceDim, self.solution.subfields.velocity.ncomponents)
    self.assertEqual(Field.VECTOR, self.solution.subfields.velocity.vectorFieldType)
    self.assertEqual(lengthScale/timeScale, self.solution.subfields.velocity.scale)

    # Fault Lagrange multiplier
    from pylith.topology.Field import Field
    self.assertEqual(spaceDim, self.solution.subfields.lagrangeFault.ncomponents)
    self.assertEqual(Field.VECTOR, self.solution.subfields.lagrangeFault.vectorFieldType)
    self.assertEqual(pressureScale, self.solution.subfields.lagrangeFault.scale)
    return


# End of file 
