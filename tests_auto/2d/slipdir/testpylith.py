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

  from TestFaultX import TestFaultX
  suite.addTest(unittest.makeSuite(TestFaultX))

  from TestFaultY import TestFaultY
  suite.addTest(unittest.makeSuite(TestFaultY))

  from TestFaultXYP import TestFaultXYP
  suite.addTest(unittest.makeSuite(TestFaultXYP))

  from TestFaultXYN import TestFaultXYN
  suite.addTest(unittest.makeSuite(TestFaultXYN))

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
