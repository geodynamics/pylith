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

## @file pyre/meshio/OutputFaultKin.py
##
## @brief Python object for managing output of finite-element
## information for faults with kinematic ruptures.
##
## Factory: output_manager

from OutputManager import OutputManager

# OutputFaultKin class
class OutputFaultKin(OutputManager):
  """
  Python object for managing output of finite-element information for
  faults with kinematic ruptures.

  Factory: output_manager
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(OutputManager.Inventory):
    """
    Python object for managing OutputFaultKin facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing OutputFaultKin facilities and properties.
    ##
    ## \b Properties
    ## @li \b vertex_info_fields Names of vertex info fields to output.
    ## @li \b vertex_data_fields Names of vertex data fields to output.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    vertexInfoFields = pyre.inventory.list("vertex_info_fields",
                                           default=["normal-dir",
                                                    "final-slip",
                                                    "slip-time"])
    vertexInfoFields.meta['tip'] = "Names of vertex info fields to output."

    vertexDataFields = pyre.inventory.list("vertex_data_fields", 
                                           default=["slip","traction-change"])
    vertexDataFields.meta['tip'] = "Names of vertex data fields to output."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="outputfaultkin"):
    """
    Constructor.
    """
    OutputManager.__init__(self, name)
    return

    
  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    OutputManager._configure(self)
    self.vertexInfoFields = self.inventory.vertexInfoFields
    self.vertexDataFields = self.inventory.vertexDataFields
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
  """
  Factory associated with OutputFaultKin.
  """
  return OutputFaultKin()


# End of file 
