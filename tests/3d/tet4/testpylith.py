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

  from TestAxialElasticIsotropic import TestAxialElasticIsotropic
  suite.addTest(unittest.makeSuite(TestAxialElasticIsotropic))

  from TestShearElasticIsotropic import TestShearElasticIsotropic
  suite.addTest(unittest.makeSuite(TestShearElasticIsotropic))

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
