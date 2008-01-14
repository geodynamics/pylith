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

## @file pyre/meshio/DataWriter.py
##
## @brief Python abstract base class for writing finite-element data.
##
## Factory: output_data_writer

from pyre.components.Component import Component

# DataWriter class
class DataWriter(Component):
  """
  Python abstract base class for writing finite-element data.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing DataWriter facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing DataWriter facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="datawriter"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="datawriter")
    self.cppHandle = None
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    return


  def initialize(self):
    """
    Initialize writer.
    """
    self._createCppHandle()

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
    Create handle to corresponding C++ object.
    """
    raise NotImplementedError("Please implement _createCppHandle() in " \
                              "derived class.")
  
  
# End of file
