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

## @file pylith/perf/Mesh.py
##
## @brief Python memory model for Mesh.

from Memory import Memory

cellTypes = {(0,1): 'point',
             (1,2): 'line2',
             (1,3): 'line3',
             (2,3): 'tri3',
             (2,6): 'tri6',
             (2,4): 'quad4',
             (2,9): 'quad9',
             (3,4): 'tet4',
             (3,8): 'hex8',
             (3,10): 'tet10',
             (3,27): 'hex27'
             }


class Mesh(Memory):
  """
  Mesh object for holding mesh memory and performance information.
  """

  def __init__(self, dimension=0, maxConeSize=0, numVertices=0, numCells=0):
    """
    Constructor.
    """
    self.dimension = dimension
    self.coneSize = maxConeSize
    self.nvertices = numVertices
    self.ncells = numCells
    self.cellType = None
    self.initialize()
    return


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
      if self.coneSize > 0:
        self.cellType = cellTypes[(self.dimension,self.coneSize)]
    except:
      raise ValueError("Unknown cell type '%s' for dim %d and cone size %d." %\
                         (self.cellType,self.dimension,self.coneSize))
    return


  def tabulate(self, memDict):
    """
    Tabulate memory use.
    """
    memDict['Creation'] = self.sizeInt * (2 * (self.coneSize*self.ncells + self.nvertices + self.ncells) + self.coneSize*self.ncells)
    memDict['Stratification'] = 2 * self.sizeArrow * (self.nvertices + self.ncells)
    # Here we have data + atlas (could use uniform) + bc (could use Section)
    memDict['Coordinates'] = (self.sizeDouble * self.dimension * self.nvertices) + (2 * self.sizeInt * self.nvertices) + (2 * self.sizeInt * self.nvertices)
    memDict['Overlap'] = 0 # Don't know overlap
    memDict['RealSections'] = 0 # Real sections should be elsewhere
    return

if __name__ == '__main__':
  d = {}
  Mesh(2, 3, 10, 25).tabulate(d)
  print 'Memory:',d
