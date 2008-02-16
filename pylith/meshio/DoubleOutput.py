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

## @file pylith/meshio/DoubleOutput.py
##
## @brief Python container with two output managers.
##
## Factory: object_bin

from pylith.utils.ObjectBin import ObjectBin

# DoubleOutput class
class DoubleOutput(ObjectBin):
  """
  Python container with two output managers.

  Factory: object_bin
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(ObjectBin.Inventory):
    """
    Python object for managing DoubleOutput facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing DoubleOutput facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b domain Output manager for domain.
    ## @li \b subdomain Output manager for subdomain.

    import pyre.inventory

    from OutputSoln import OutputSoln
    domain = pyre.inventory.facility("domain", family="output_manager",
                                     factory=OutputSoln)
    domain.meta['tip'] = "Output manager for domain."

    from OutputSolnSubset import OutputSolnSubset
    subdomain = pyre.inventory.facility("subdomain", family="output_manager",
                                        factory=OutputSolnSubset)
    subdomain.meta['tip'] = "Output manager for subdomain."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="doubleoutput"):
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
    self.bin = [self.inventory.domain,
                self.inventory.subdomain]
    return

  
# FACTORIES ////////////////////////////////////////////////////////////

def object_bin():
  """
  Factory associated with DoubleOutput.
  """
  return DoubleOutput()


# End of file 
