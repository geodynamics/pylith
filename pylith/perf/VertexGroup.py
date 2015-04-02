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

## @file pylith/perf/VertexGroup.py
##
## @brief Python memory model for vertex groups.

from Memory import Memory

class VertexGroup(Memory):
  """
  Mesh object for holding vertex group memory and performance information.
  """
  def __init__(self, label = '', numVertices = 0, numMeshVertices = 0):
    """
    Constructor.
    """
    self.label         = label
    self.nvertices     = numVertices
    self.nMeshVertices = numMeshVertices
    return

  def tabulate(self, memDict):
    """
    Tabulate memory use.
    """
    # Here we have data + atlas (could use uniform) + bc (use Section)
    if not self.label in memDict:
      memDict[self.label] = 0
    memDict[self.label] += (self.sizeInt * self.nvertices) + (2 * self.sizeInt * self.nMeshVertices) + (2 * self.sizeInt * self.nMeshVertices)
    return

if __name__ == '__main__':
  d = {}
  VertexGroup('rock', 35).tabulate(d)
  print 'Memory:',d


# End of file

