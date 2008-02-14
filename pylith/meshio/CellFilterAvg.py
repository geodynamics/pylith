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

## @file pyre/meshio/CellFilterAvg.py
##
## @brief Python class for averageing cell fields over each cell's
## quadrature points when writing finite-element data.
##
## Factory: output_cell_filter

from CellFilter import CellFilter

# CellFilterAvg class
class CellFilterAvg(CellFilter):
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
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _createCppHandle(self):
    """
    Create handle to C++ object.
    """
    if None == self.cppHandle:
      import meshio as bindings
      self.cppHandle = bindings.CellFilterAvg()
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def output_cell_filter():
  """
  Factory associated with CellFilter.
  """
  return CellFilterAvg()


# End of file 
