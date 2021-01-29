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
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file tests/pytests/problems/TestProgressMonitor.py

## @brief Unit testing of ProgressMonitor object.

import unittest
from pylith.problems.ProgressMonitor import ProgressMonitor

from pythia.pyre.units.time import year

# ----------------------------------------------------------------------
class TestProgressMonitor(unittest.TestCase):
  """
  Unit testing of ProgressMonitor object.
  """

  def setUp(self):
    self.monitor = ProgressMonitor()
    self.monitor._configure()
    return
  

  def test_constructor(self):
    """
    Test constructor.
    """
    monitor = ProgressMonitor()
    monitor._configure()
    return


  def test_openclose(self):
    """
    Test open() and close().
    """
    self.monitor.open()
    self.monitor.close()
    return


  def test_update(self):
    """
    Test update().
    """
    self.monitor.open()

    self.monitor.update(1.0*year, 0.0*year, 10.0*year)
    self.monitor.update(1.5*year, 0.0*year, 10.0*year)
    self.monitor.update(2.0*year, 0.0*year, 10.0*year)
    self.monitor.update(4.0*year, 0.0*year, 10.0*year)
    self.monitor.update(5.0*year, 0.0*year, 10.0*year)
    self.monitor.update(6.0*year, 0.0*year, 10.0*year)
    self.monitor.update(8.0*year, 0.0*year, 10.0*year)
    self.monitor.update(9.0*year, 0.0*year, 10.0*year)

    self.monitor.close()

    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.problems.ProgressMonitor import progress_monitor
    m = progress_monitor()
    return


# End of file 
