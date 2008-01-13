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

## @file pyre/meshio/VertexFilter.py
##
## @brief Python abstract base class for filtering vertex fields when
## writing finite-element data.
##
## Factory: output_vertex_filter

from pyre.components.Component import Component

# VertexFilter class
class VertexFilter(Component):
  """
  Python abstract base class for filtering cell fields when writing
  finite-element data.

  Factory: output_vertex_filter
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing VertexFilter facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing VertexFilter facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="vertexfilter"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="vertexfilter")
    self.cppHandle = None
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    return


  def initialize(self):
    """
    Initialize output manager.
    """
    if None == self.cppHandle:
      import pylith.meshio.meshio as bindings
      self.cppHandle = bindings.VertexFilter()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_vertex_filter():
  """
  Factory associated with VertexFilter.
  """
  return VertexFilter()


# End of file 
