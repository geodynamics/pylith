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
## @brief Python container with one output type.
##
## Factory: object_bin

from pylith.utils.ObjectBin import ObjectBin

# SingleOutput class
class SingleOutput(ObjectBin):
  """
  Python container with one output type.

  Factory: object_bin
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(ObjectBin.Inventory):
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

    from SolutionIOVTK import SolutionIOVTK
    output = pyre.inventory.facility("output", family="solution_io",
                                     factory=SolutionIOVTK)
    output.meta['tip'] = "Output manager."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="singleoutput"):
    """
    Constructor.
    """
    ObjectBin.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set attributes from inventory.
    """
    ObjectBin._configure(self)
    self.bin = [self.inventory.output]
    return

  
# FACTORIES ////////////////////////////////////////////////////////////

def object_bin():
  """
  Factory associated with SingleOutput.
  """
  return SingleOutput()


# End of file 
