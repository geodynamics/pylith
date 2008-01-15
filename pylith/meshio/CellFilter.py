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

from pyre.components.Component import Component

# CellFilter class
class CellFilter(Component):
  """
  Python abstract base class for filtering cell fields when writing
  finite-element data.

  Factory: output_cell_filter
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing CellFilter facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing CellFilter facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="cellfilter"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="cellfilter")
    self.cppHandle = None
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
    self._createCppHandle()

    if None != self.cppHandle and quadrature != None:
      # Only set quadrature if filter is specified and quadrature is
      # provided.
      assert(None != self.cppHandle)
      self.cppHandle.quadrature = quadrature.cppHandle
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    return


  def _createCppHandle(self):
    """
    Create handle to C++ object.
    """
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def output_cell_filter():
  """
  Factory associated with CellFilter.
  """
  return CellFilter()


# End of file 
