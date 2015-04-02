#!/usr/bin/env nemesis
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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

import unittest

def suite():
  """
  Create test suite.
  """
  suite = unittest.TestSuite()

  from TestAxialDisp import TestAxialDisp
  suite.addTest(unittest.makeSuite(TestAxialDisp))

  from TestShearDisp import TestShearDisp
  suite.addTest(unittest.makeSuite(TestShearDisp))

  from TestDislocation import TestDislocation
  suite.addTest(unittest.makeSuite(TestDislocation))

  from TestDislocation import TestDislocation2
  suite.addTest(unittest.makeSuite(TestDislocation2))

  # Not complete
  ##from TestDislocationTwoFaults import TestDislocation
  ##suite.addTest(unittest.makeSuite(TestDislocation))

  from TestLgDeformRigidBody import TestRigidBody
  suite.addTest(unittest.makeSuite(TestRigidBody))

  from TestLgDeformTraction import TestTraction
  suite.addTest(unittest.makeSuite(TestTraction))

  from TestFrictionCompression import TestFrictionCompression
  suite.addTest(unittest.makeSuite(TestFrictionCompression))
  
  from TestFrictionOpening import TestFrictionOpening
  suite.addTest(unittest.makeSuite(TestFrictionOpening))

  from TestFrictionShearStick import TestFrictionShearStick
  suite.addTest(unittest.makeSuite(TestFrictionShearStick))

  from TestFrictionShearSliding import TestFrictionShearSliding
  suite.addTest(unittest.makeSuite(TestFrictionShearSliding))

  from TestSlipWeakeningCompression import TestSlipWeakeningCompression
  suite.addTest(unittest.makeSuite(TestSlipWeakeningCompression))
  
  from TestSlipWeakeningOpening import TestSlipWeakeningOpening
  suite.addTest(unittest.makeSuite(TestSlipWeakeningOpening))

  from TestSlipWeakeningShearStick import TestSlipWeakeningShearStick
  suite.addTest(unittest.makeSuite(TestSlipWeakeningShearStick))

  from TestSlipWeakeningShearSliding import TestSlipWeakeningShearSliding
  suite.addTest(unittest.makeSuite(TestSlipWeakeningShearSliding))

  return suite


def main():
  """
  Run test suite.
  """
  unittest.TextTestRunner(verbosity=2).run(suite())
  return


# ----------------------------------------------------------------------
if __name__ == '__main__':
  main()

  
# End of file 
