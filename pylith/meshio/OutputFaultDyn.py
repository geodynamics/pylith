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

## @file pyre/meshio/OutputFaultDyn.py
##
## @brief Python object for managing output of finite-element
## information for faults with dynamic ruptures.
##
## Factory: output_manager

from OutputManager import OutputManager

# OutputFaultDyn class
class OutputFaultDyn(OutputManager):
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

  cellInfoFields = pyre.inventory.list("cell_info_fields", default=[])
  cellInfoFields.meta['tip'] = "Names of cell info fields to output."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="outputfaultdyn"):
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
    self.cellInfoFields   = self.inventory.cellInfoFields
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
  """
  Factory associated with OutputManager.
  """
  return OutputFaultDyn()


# End of file 
