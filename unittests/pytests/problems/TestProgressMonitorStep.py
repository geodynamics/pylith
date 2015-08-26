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

## @file unittests/pytests/problems/TestProgressMonitorStep.py

## @brief Unit testing of ProgressMonitorStep object.

import unittest
from pylith.problems.ProgressMonitorStep import ProgressMonitorStep

from pyre.units.time import year

# ----------------------------------------------------------------------
class TestProgressMonitorStep(unittest.TestCase):
  """
  Unit testing of ProgressMonitorStep object.
  """

  def setUp(self):
    self.monitor = ProgressMonitorStep()
    self.monitor._configure()
    self.monitor.filename = "data/progress_step.txt"
    return
  

  def test_constructor(self):
    """
    Test constructor.
    """
    monitor = ProgressMonitorStep()
    monitor._configure()
    return


  def test_openclose(self):
    """
    Test open() and close().
    """
    import os
    if os.path.exists(self.monitor.filename):
        os.remove(self.monitor.filename)
    self.monitor.open()
    self.monitor.close()

    self.assertTrue(os.path.isfile(self.monitor.filename))

    return


  def test_update(self):
    """
    Test update().
    """
    import os
    self.monitor.open()

    nlines = 1 # header
    self.monitor.update(1, 1, 100); nlines += 1
    self.monitor.update(2, 1, 100); nlines += 0
    self.monitor.update(3, 1, 100); nlines += 0
    self.monitor.update(6, 1, 100); nlines += 1
    self.monitor.update(10, 1, 100); nlines += 1
    self.monitor.update(12, 1, 100); nlines += 0
    self.monitor.update(20, 1, 100); nlines += 1
    self.monitor.update(40, 1, 100); nlines += 1
    self.monitor.update(50, 1, 100); nlines += 1
    self.monitor.update(60, 1, 100); nlines += 1
    self.monitor.update(80, 1, 100); nlines += 1
    self.monitor.update(90, 1, 100); nlines += 1

    self.monitor.close()

    self.assertTrue(os.path.isfile(self.monitor.filename))
    fin = open(self.monitor.filename, "r")
    lines = fin.readlines()
    fin.close()
    self.assertEqual(nlines, len(lines))
    
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.problems.ProgressMonitorStep import progress_monitor
    m = progress_monitor()
    return


# End of file 
