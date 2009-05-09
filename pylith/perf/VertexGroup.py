#!/usr/bin/env python

from Memory import Memory

class VertexGroup(Memory):
  """
  Mesh object for holding vertex group memory and performance information.
  """
  def __init__(self, label = '', numVertices = 0):
    """
    Constructor.
    """
    self.label     = label
    self.nvertices = numVertices
    return

  def tabulateNew(self, memDict):
    """
    Tabulate memory use.
    """
    memDict[self.label] = self.sizeInt * (2 * self.nvertices + self.nvertices)
    return

  def tabulate(self):
    """
    Tabulate memory use.
    """
    return self.sizeInt * ( 2 * self.nvertices + self.nvertices)

if __name__ == '__main__':
  print 'Memory:',VertexGroup('rock', 35).tabulate()
