#!/usr/bin/env python

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
  print 'Memory:',VertexGroup('rock', 35).tabulate()
