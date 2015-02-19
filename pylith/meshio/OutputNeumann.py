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

## @file pyre/meshio/OutputNeumann.py
##
## @brief Python object for managing output of finite-element
## information for Neumann boundary conditions.
##
## Factory: output_manager

from OutputManager import OutputManager

# OutputNeumann class
class OutputNeumann(OutputManager):
  """
  Python object for managing output of finite-element information for
  Neumann boundary conditions.

  Factory: output_manager
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(OutputManager.Inventory):
    """
    Python object for managing OutputNeumann facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing OutputNeumann facilities and properties.
    ##
    ## \b Properties
    ## @li \b cell_info_fields Names of cell info fields to output.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    cellInfoFields = pyre.inventory.list("cell_info_fields",
                                         default=[])
    cellInfoFields.meta['tip'] = "Names of cell info fields to output."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="outputneumann"):
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
    self.vertexInfoFields = []
    self.vertexDataFields = []
    self.cellInfoFields = self.inventory.cellInfoFields
    self.cellDataFields = []
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
  """
  Factory associated with OutputManager.
  """
  return OutputNeumann()


# End of file 
