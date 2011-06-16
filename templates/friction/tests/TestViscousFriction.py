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

# We cannot test the low-level functionality of the ViscousFriction
# object because it is not exposed to Python. You should really setup
# C++ unit tests using CppUnit as is done for PyLith in addition to
# the simple Python unit tests here.

import unittest


class TestViscousFriction(unittest.TestCase):
  """
  Unit testing of ViscousFriction object.
  """

  def setUp(self):
    """
    Setup test subject.
    """
    from pylith.friction.contrib.ViscousFriction import ViscousFriction
    self.model = ViscousFriction()
    return
  

  def test_label(self):
    """
    Test constructor.
    """
    label = "viscous friction"
    self.model.label(label)
    self.assertEqual(label, self.model.label())
    return


  def test_timeStep(self):
    """
    Test constructor.
    """
    dt = 2.4
    self.model.timeStep(dt)
    self.assertEqual(dt, self.model.timeStep())
    return


  def testHasProperty(self):
    self.failUnless(self.model.hasProperty("static_coefficient"))
    self.failUnless(self.model.hasProperty("reference_slip_rate"))
    self.failUnless(self.model.hasProperty("cohesion"))
    return


  def testHasStateVar(self):
    self.failUnless(self.model.hasStateVar("slip_rate"))
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.friction.contrib.ViscousFriction import friction_model
    f = friction_model()
    return


# End of file 
