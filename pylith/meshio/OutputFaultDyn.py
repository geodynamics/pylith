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

## @file pyre/meshio/OutputFaultDyn.py
##
## @brief Python object for managing output of finite-element
## information for faults with dynamic ruptures.
##
## Factory: output_manager

from OutputManagerSubMesh import OutputManagerSubMesh

# OutputFaultDyn class
class OutputFaultDyn(OutputManagerSubMesh):
  """
  Python object for managing output of finite-element information for
  faults with dynamic ruptures.

  Inventory

  @class Inventory
  Python object for managing OutputFaultDyn facilities and properties.
  
  \b Properties
  @li \b cell_info_fields Names of cell info fields to output.
  @li \b cell_data_fields Names of cell data fields to output.
  
  \b Facilities
  @li None

  Factory: output_manager
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  cellInfoFields = pyre.inventory.list("cell_info_fields",
                                         default=["normal_dir"])
  cellInfoFields.meta['tip'] = "Names of cell info fields to output."

  cellDataFields = pyre.inventory.list("cell_data_fields", 
                                       default=["slip",
                                                "traction"])
  cellDataFields.meta['tip'] = "Names of cell data fields to output."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="outputfaultdyn"):
    """
    Constructor.
    """
    OutputManagerSubMesh.__init__(self, name)
    return

    
  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    OutputManagerSubMesh._configure(self)
    self.cellInfoFields = self.inventory.cellInfoFields
    self.cellDataFields = self.inventory.cellDataFields
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
  """
  Factory associated with OutputManager.
  """
  return OutputFaultDyn()


# End of file 
