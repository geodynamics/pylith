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

## @file pylith/meshio/SingleOutput.py
##
## @brief Python container with one output manager.

from pyre.components.Component import Component

# SingleOutput class
class SingleOutput(Component):
  """
  Python container with one output manager.

  Factory: object_bin
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing SingleOutput facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing SingleOutput facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b output Output manager

    import pyre.inventory

    from OutputSoln import OutputSoln
    output = pyre.inventory.facility("output", family="output_manager",
                                     factory=OutputSoln)
    output.meta['tip'] = "Output manager."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="singleoutput"):
    """
    Constructor.
    """
    Component.__init__(self, name)
    return


# End of file 
