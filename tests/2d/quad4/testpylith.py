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

__requires__ = "PyLith"

import unittest

def suite():
  """
  Create test suite.
  """
  suite = unittest.TestSuite()

  from TestAxialShearDisp import TestAxialShearDisp
  suite.addTest(unittest.makeSuite(TestAxialShearDisp))

  #from TestAxialShearTract import TestAxislShearTract
  #suite.addTest(unittest.makeSuite(TestAxialShearTract))

  #from TestDislocation import TestDislocation
  #suite.addTest(unittest.makeSuite(TestDislocation))

  #from TestDislocation2 import TestDislocation2
  #suite.addTest(unittest.makeSuite(TestDislocation2))

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
