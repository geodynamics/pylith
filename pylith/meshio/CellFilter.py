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

## @file pyre/meshio/CellFilter.py
##
## @brief Python abstract base class for filtering cell fields when
## writing finite-element data.
##
## Factory: output_cell_filter

from pylith.utils.PetscComponent import PetscComponent

# CellFilter class
class CellFilter(PetscComponent):
  """
  Python abstract base class for filtering cell fields when writing
  finite-element data.

  Factory: output_cell_filter
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  from pylith.perf.MemoryLogger import MemoryLogger
  perfLogger = pyre.inventory.facility("perf_logger", family="perf_logger", factory=MemoryLogger)
  perfLogger.meta['tip'] = "Performance and memory logging."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="cellfilter"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="cellfilter")
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    return


  def initialize(self, quadrature):
    """
    Initialize output manager.
    """
    return


  def finalize(self):
    """
    Cleanup after running problem.
    """
    self._modelMemoryUse()
    return


  # PRIVATE METHODS ///////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)
    self.perfLogger = self.inventory.perfLogger
    return


  def _modelMemoryUse(self):
    """
    Model memory allocation.
    """
    raise NotImplementedError("Please implement _modelMemoryUse() in derived class.")
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_cell_filter():
  """
  Factory associated with CellFilter.
  """
  # Abstract object (so return None).
  return None


# End of file 
