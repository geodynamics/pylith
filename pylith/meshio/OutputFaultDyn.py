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
  @li \b vertex_info_fields Names of vertex info fields to output.
  @li \b vertex_data_fields Names of vertex data fields to output.
  
  \b Facilities
  @li None

  Factory: output_manager
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  vertexInfoFields = pyre.inventory.list("vertex_info_fields",
                                         default=["normal_dir"])
  vertexInfoFields.meta['tip'] = "Names of vertex info fields to output."

  vertexDataFields = pyre.inventory.list("vertex_data_fields", 
                                         default=["slip",
                                                  "traction"])
  vertexDataFields.meta['tip'] = "Names of vertex data fields to output."

  cellInfoFields = pyre.inventory.list("cell_info_fields",
                                       default=["distribution"])
  cellInfoFields.meta['tip'] = "Names of cell info fields to output."


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
    self.vertexInfoFields = self.inventory.vertexInfoFields
    self.vertexDataFields = self.inventory.vertexDataFields
    self.cellInfoFields   = self.inventory.cellInfoFields
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
  """
  Factory associated with OutputManager.
  """
  return OutputFaultDyn()


# End of file 
