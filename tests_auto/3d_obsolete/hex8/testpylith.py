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

  from TestAxialElasticIsotropic import TestAxialElasticIsotropic
  suite.addTest(unittest.makeSuite(TestAxialElasticIsotropic))

  from TestShearElasticIsotropic import TestShearElasticIsotropic
  suite.addTest(unittest.makeSuite(TestShearElasticIsotropic))

  from TestShearPlaneStrain import TestShearPlaneStrain
  suite.addTest(unittest.makeSuite(TestShearPlaneStrain))

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
