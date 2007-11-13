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

  from TestAxial import TestAxial
  suite.addTest(unittest.makeSuite(TestAxial))

  from TestDislocationStatic import TestDislocationStatic
  suite.addTest(unittest.makeSuite(TestDislocationStatic))

  from TestDislocationDyn import TestDislocationDyn
  suite.addTest(unittest.makeSuite(TestDislocationDyn))

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
