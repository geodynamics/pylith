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

## @file pyre/meshio/OutputDirichlet.py
##
## @brief Python object for managing output of finite-element
## information for Dirichlet boundary conditions.
##
## Factory: output_manager

from OutputManager import OutputManager

# OutputDirichlet class
class OutputDirichlet(OutputManager):
  """
  Python object for managing output of finite-element information for
  Dirichlet boundary conditions.

  Factory: output_manager
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(OutputManager.Inventory):
    """
    Python object for managing OutputDirichlet facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing OutputDirichlet facilities and properties.
    ##
    ## \b Properties
    ## @li \b vertex_info_fields Names of vertex info fields to output.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    vertexInfoFields = pyre.inventory.list("vertex_info_fields",
                                           default=[])
    vertexInfoFields.meta['tip'] = "Names of vertex info fields to output."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="outputdirichlet"):
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
    self.vertexDataFields = []
    self.cellInfoFields = []
    self.cellDataFields = []
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
  """
  Factory associated with OutputManager.
  """
  return OutputDirichlet()


# End of file 
