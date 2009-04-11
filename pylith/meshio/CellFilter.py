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


# FACTORIES ////////////////////////////////////////////////////////////

def output_cell_filter():
  """
  Factory associated with CellFilter.
  """
  return CellFilter()


# End of file 
