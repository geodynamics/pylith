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

  from TestAxialDisp import TestAxialDisp
  suite.addTest(unittest.makeSuite(TestAxialDisp))

  from TestShearDisp import TestShearDisp
  suite.addTest(unittest.makeSuite(TestShearDisp))

  from TestDislocation import TestDislocation
  suite.addTest(unittest.makeSuite(TestDislocation))

  #from TestDislocation2 import TestDislocation2
  #suite.addTest(unittest.makeSuite(TestDislocation2))

  from TestLgDeformRigidBody import TestRigidBody
  suite.addTest(unittest.makeSuite(TestRigidBody))

  from TestLgDeformTraction import TestTraction
  suite.addTest(unittest.makeSuite(TestTraction))

  from TestFrictionCompression import TestFrictionCompression
  suite.addTest(unittest.makeSuite(TestFrictionCompression))
  
  from TestFrictionShearStick import TestFrictionShearStick
  suite.addTest(unittest.makeSuite(TestFrictionShearStick))

  from TestFrictionShearSliding import TestFrictionShearSliding
  suite.addTest(unittest.makeSuite(TestFrictionShearSliding))

  from TestFrictionOpening import TestFrictionOpening
  suite.addTest(unittest.makeSuite(TestFrictionOpening))

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
