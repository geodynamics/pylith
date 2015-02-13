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

## @file pylith/perf/Field.py
##
## @brief Python memory model for field.

from Memory import Memory

class Field(Memory):
  """
  Mesh object for holding field memory and performance information.
  """
  def __init__(self, label = '', size = 0, chartSize = 0):
    """
    Constructor.
    """
    self.label     = label
    self.size      = size
    self.chartSize = chartSize
    return


  def tabulate(self, memDict):
    """
    Tabulate memory use.
    """
    # Here we have data + atlas + bc
    if not self.label in memDict:
      memDict[self.label] = 0
    memDict[self.label] += \
        (self.sizeDouble * self.size) + \
        (2 * self.sizeInt * self.chartSize) + \
        (2 * self.sizeInt * self.chartSize)
    return


if __name__ == '__main__':
  print 'Memory:',Material('rock', 35).tabulate()


# End of file
