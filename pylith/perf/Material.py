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

## @file pylith/perf/Material.py
##
## @brief Python memory model for material.

from Memory import Memory

class Material(Memory):
  """
  Mesh object for holding material memory and performance information.
  """
  def __init__(self, label = '', numCells = 0):
    """
    Constructor.
    """
    self.label  = label
    self.ncells = numCells
    return

  def tabulate(self, memDict):
    """
    Tabulate memory use.
    """
    if not self.label in memDict:
      memDict[self.label] = 0
    memDict[self.label] += self.sizeArrow * self.ncells
    return

if __name__ == '__main__':
  print 'Memory:',Material('rock', 35).tabulate()


# End of file
