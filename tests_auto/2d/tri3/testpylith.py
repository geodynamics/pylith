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
# Copyright (c) 2010-2014 University of California, Davis
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

  from TestShearDispNoSlip import TestShearDispNoSlip
  suite.addTest(unittest.makeSuite(TestShearDispNoSlip))

  from TestShearDispFriction import TestShearDispFriction
  suite.addTest(unittest.makeSuite(TestShearDispFriction))

  from TestShearDispNoSlipRefine import TestShearDispNoSlipRefine
  suite.addTest(unittest.makeSuite(TestShearDispNoSlipRefine))

  from TestSlipOneFault import TestSlipOneFault
  suite.addTest(unittest.makeSuite(TestSlipOneFault))

  from TestSlipTwoFaults import TestSlipTwoFaults
  suite.addTest(unittest.makeSuite(TestSlipTwoFaults))

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
