#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

## @file pylith/perf/Jacobian.py
##
## @brief Python memory model for Jacobian sparse matrix.

from Memory import Memory

class Jacobian(Memory):
  """
  Mesh object for holding matrix memory and performance information.
  """
  def __init__(self, label = ''):
    """
    Constructor.
    """
    self.label = label
    return


  def tabulate(self, memDict):
    """
    Tabulate memory use.
    """
    # Here we have preallocation sieve and Mat
    if not self.label in memDict:
      memDict[self.label] = 0
    memDict[self.label] += 0
    return


# End of file
