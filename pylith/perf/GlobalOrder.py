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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/perf/GlobalOrder.py
##
## @brief Python memory model for global order.

from Memory import Memory

class GlobalOrder(Memory):
  """
  Mesh object for holding global order memory and performance information.
  """
  def __init__(self, label = '', chartSize = 0):
    """
    Constructor.
    """
    self.label     = label
    self.chartSize = chartSize
    return


  def tabulate(self, memDict):
    """
    Tabulate memory use.
    """
    # Here we have a map<int --> (int,int)> + set<int>
    if not self.label in memDict:
      memDict[self.label] = 0
    memDict[self.label] += self.chartSize*(3 * self.sizeInt + self.sizeMapEntry) + \
        self.chartSize*(self.sizeSetEntry+self.sizeInt)
    return


# End of file
