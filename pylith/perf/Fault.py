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
# Copyright (c) 2010-2011 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/perf/Fault.py
##
## @brief Python memory model for Fault.

from Memory import Memory

class Fault(Memory):
  """
  Fault object for holding mesh memory and performance information.
  """
  def __init__(self):
    """
    Constructor.
    """
    return


  def tabulate(self, memDict):
    """
    Tabulate memory use.
    """
    memDict['Creation'] = 0
    memDict['Stratification'] = 0
    memDict['Coordinates'] = 0
    return

