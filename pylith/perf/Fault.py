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

## @file pylith/perf/Fault.py
##
## @brief Python memory model for Fault.

from Memory import Memory

class Fault(Memory):
  """
  Fault object for holding mesh memory and performance information.
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


  def initialize(self):
    """
    Initialize application.
    """
    from Mesh import cellTypes
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
    dim = self.dimension
    ncells = self.ncells
    nvertices = self.nvertices
    coneSize = self.coneSize

    memDict['Creation'] = self.sizeInt * (2 * (coneSize*ncells + nvertices + ncells) + coneSize*ncells)
    memDict['Stratification'] = 2 * self.sizeArrow * (nvertices + ncells)
    memDict['Coordinates'] = (self.sizeDouble * dim * nvertices) + (2 * self.sizeInt * nvertices) + (2 * self.sizeInt * nvertices)

    return

