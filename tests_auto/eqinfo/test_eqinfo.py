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

  from TestEqInfoLine import TestEqInfoLine
  suite.addTest(unittest.makeSuite(TestEqInfoLine))

  from TestEqInfoTri3 import TestEqInfoTri3
  suite.addTest(unittest.makeSuite(TestEqInfoTri3))

  from TestEqInfoQuad4 import TestEqInfoQuad4
  suite.addTest(unittest.makeSuite(TestEqInfoQuad4))

  return suite


def main():
  """
  Run test suite.
  """
  import sys
  sys.path.append(".")

  unittest.TextTestRunner(verbosity=2).run(suite())
  return


# ----------------------------------------------------------------------
if __name__ == '__main__':
  main()

  
# End of file 
