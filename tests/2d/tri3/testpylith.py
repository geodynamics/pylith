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

  from TestAxialPlaneStrain import TestAxialPlaneStrain
  suite.addTest(unittest.makeSuite(TestAxialPlaneStrain))

  from TestShearPlaneStrain import TestShearPlaneStrain
  suite.addTest(unittest.makeSuite(TestShearPlaneStrain))

  from TestDislocation import TestDislocation
  suite.addTest(unittest.makeSuite(TestDislocation))

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
