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

## @file pyre/meshio/CellFilterAvg.py
##
## @brief Python class for averaging cell fields over each cell's
## quadrature points when writing finite-element data.
##
## Factory: output_cell_filter

from CellFilter import CellFilter
from meshio import CellFilterAvg as ModuleCellFilterAvg

# CellFilterAvg class
class CellFilterAvg(CellFilter, ModuleCellFilterAvg):
  """
  Python class for average cell fields over each cell's quadrature
  points when writing finite-element data.

  Factory: output_cell_filter
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="cellfilteravg"):
    """
    Constructor.
    """
    CellFilter.__init__(self, name)
    ModuleCellFilterAvg.__init__(self)
    return


  def initialize(self, quadrature):
    """
    Initialize output manager.
    """
    self.quadrature(quadrature)
    return


  # PRIVATE METHODS ///////////////////////////////////////////////////

  def _modelMemoryUse(self):
    """
    Model memory allocation.
    """
    self.perfLogger.logField('Output', self.fieldAvg())
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_cell_filter():
  """
  Factory associated with CellFilter.
  """
  return CellFilterAvg()


# End of file 
