#!/usr/bin/env python

from Memory import Memory

class Mesh(Memory):
  cellTypes = {(1,2): 'line2',
               (2,3): 'tri3',
               (2,4): 'quad4',
               (3,4): 'tet4',
               (3,8): 'hex8'
               }

  """
  Mesh object for holding mesh memory and performance information.
  """
  def __init__(self, dimension = 0, maxConeSize = 0, numVertices = 0, numCells = 0):
    """
    Constructor.
    """
    self.dimension = dimension
    self.coneSize  = maxConeSize
    self.nvertices = numVertices
    self.ncells    = numCells
    self.cellType  = None
    self.initialize()
    return

  @classmethod
  def cellTypeInfo(cls, cellType):
    for k,cT in cls.cellTypes.iteritems():
      if cT == cellType:
        return k
    raise ValueError("Unknown cell type '%s'." % cellType)

  def initialize(self):
    """
    Initialize application.
    """
    try:
      self.cellType = self.cellTypes[(self.dimension,self.coneSize)]
    except:
      raise ValueError("Unknown cell type '%s' for dim %d and cone size %d." % (self.cellType,self.dimension,self.coneSize))
    return

  def tabulate(self, memDict):
    """
    Tabulate memory use.
    """
    memDict['Creation']       = self.sizeInt * (2 * (self.coneSize*self.ncells + self.nvertices + self.ncells) + self.coneSize*self.ncells)
    memDict['Stratification'] = 2 * self.sizeArrow * (self.nvertices + self.ncells)
    # Here we have data + atlas (could use uniform) + bc (could use Section)
    memDict['Coordinates']    = (self.sizeDouble * self.dimension * self.nvertices) + (2 * self.sizeInt * self.nvertices) + (2 * self.sizeInt * self.nvertices)
    return

if __name__ == '__main__':
  print 'Memory:',Mesh(2, 3, 10, 25).tabulate()
