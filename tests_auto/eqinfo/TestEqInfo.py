#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2012 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file tests_auto/eqinfo/TestEqInfo.py
##
## @brief Generic tests for pylith_eqinfo.

import unittest
import numpy

class TestEqInfo(unittest.TestCase):
  """
  Generic tests for pylith_eqinfo.
  """

  def _check(self, statsE, stats):
    """
    Check earthquake stats.
    """
    tol = 1.0e-6
    attrs = ["timestamp",
             "ruparea",
             "potency",
             "moment",
             "avgslip",
             "mommag",
             ]

    statsE.avgslip = statsE.potency / (statsE.ruparea + 1.0e-30)
    statsE.mommag = 2.0/3.0*(numpy.log10(statsE.moment) - 9.05)

    for attr in attrs:
      valuesE = statsE.__getattribute__(attr)
      values = stats.__getattribute__(attr)
      msg = "Mismatch in number of snapshots for attribute '%s', %d != %d." % (attr, len(valuesE), len(values))
      self.assertEqual(len(valuesE), len(values), msg=msg)
      
      for (valueE, value) in zip(valuesE, values):
        msg = "Mismatch in value for attribute '%s', %g != %g." % (attr, valueE, value)
        if valueE != 0.0:
          self.assertAlmostEqual(valueE, value, msg=msg, delta=abs(tol*valueE))
        else:
          self.assertAlmostEqual(valueE, value, msg=msg, delta=tol)
        
    return


# End of file
