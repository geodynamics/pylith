#!/usr/bin/env python

from Memory import Memory

class Mesh(Memory):
  cellTypes = {(2,3): 'tri3',
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

  def tabulate(self):
    """
    Tabulate memory use.
    """
    memory = {'mesh': 0,
              'stratification': 0,
              'coordinates': 0,
              'materials': 0}
    ncells    = self.ncells
    nvertices = self.nvertices
    coneSize  = self.coneSize
    dimension = self.dimension

    # mesh
    nbytes = self.sizeInt * (2 * (coneSize*ncells + nvertices + ncells) + coneSize*ncells)
    memory['mesh'] = nbytes

    # stratification
    nbytes = 2 * self.sizeArrow * (nvertices + ncells)
    memory['stratification'] = nbytes

    # coordinates
    nbytes = self.sizeDouble * dimension * nvertices
    memory['coordinates'] = nbytes

    # materials
    nbytes = 2 * self.sizeArrow * ncells
    memory['materials'] = nbytes
    return memory

if __name__ == '__main__':
  print 'Memory:',Mesh(2, 3, 10, 25).tabulate()
