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

## @file unittests/pytests/problems/TestProgressMonitorTime.py

## @brief Unit testing of ProgressMonitorTime object.

import unittest
from pylith.problems.ProgressMonitorTime import ProgressMonitorTime

from pyre.units.time import year

# ----------------------------------------------------------------------
class TestProgressMonitorTime(unittest.TestCase):
  """
  Unit testing of ProgressMonitorTime object.
  """

  def setUp(self):
    self.monitor = ProgressMonitorTime()
    self.monitor._configure()
    self.monitor.filename = "data/progress_time.txt"
    return
  

  def test_constructor(self):
    """
    Test constructor.
    """
    monitor = ProgressMonitorTime()
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
    self.monitor.update(0.0*year, 0.0*year, 10.0*year); nlines += 1
    self.monitor.update(0.5*year, 0.0*year, 10.0*year); nlines += 1
    self.monitor.update(1.0*year, 0.0*year, 10.0*year); nlines += 1
    self.monitor.update(1.1*year, 0.0*year, 10.0*year); nlines += 0
    self.monitor.update(2.0*year, 0.0*year, 10.0*year); nlines += 1
    self.monitor.update(4.0*year, 0.0*year, 10.0*year); nlines += 1
    self.monitor.update(5.0*year, 0.0*year, 10.0*year); nlines += 1
    self.monitor.update(6.0*year, 0.0*year, 10.0*year); nlines += 1
    self.monitor.update(8.0*year, 0.0*year, 10.0*year); nlines += 1
    self.monitor.update(9.0*year, 0.0*year, 10.0*year); nlines += 1

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
    from pylith.problems.ProgressMonitorTime import progress_monitor
    m = progress_monitor()
    return


# End of file 
